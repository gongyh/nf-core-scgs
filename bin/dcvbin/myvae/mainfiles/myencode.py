from typing import Optional, IO, Union
from pathlib import Path
import vambtools as _vambtools
from torch.utils.data.dataset import TensorDataset as _TensorDataset
from torch.utils.data import DataLoader as _DataLoader
from torch.nn.functional import softmax as _softmax
from torch.optim import Adam as _Adam
from torch import Tensor
from torch import nn as _nn
from math import log as _log
from loguru import logger
import sys
from datetime import datetime
import time

__doc__ = """Encode a depths matrix and a tnf matrix to latent representation.

Creates a variational autoencoder in PyTorch and tries to represent the depths
and tnf in the latent space under gaussian noise.

Usage:
>>> vae = VAE(nsamples=6)
>>> dataloader = make_dataloader(depths, tnf, lengths)
>>> vae.trainmodel(dataloader)
>>> latent = vae.encode(dataloader) # Encode to latent representation
>>> latent.shape
(183882, 32)
"""

__cmd_doc__ = """Encode depths and TNF using a VAE to latent representation"""

import numpy as _np
import torch as _torch

# 设置日志格式
def format_log(record) -> str:
    colors = {"WARNING": "red", "INFO": "green", "DEBUG": "blue", "ERROR": "red"}
    L = colors.get(record["level"].name, "blue")
    T = "red" if record["level"].name in ("WARNING", "ERROR") else "cyan"
    message = "<red>{message}</red>" if T == "red" else "{message}"
    time = f"<{T}>{{time:YYYY-MM-DD HH:mm:ss.SSS}}</{T}>"
    level = f"<b><{L}>{{level:<7}}</{L}></b>"
    return f"{time} | {level} | {message}\n{{exception}}"
# 清空所有设置
logger.remove()
# 添加控制台输出的格式,sys.stdout为输出到屏幕;关于这些配置还需要自定义请移步官网查看相关参数说明
logger.add(sys.stdout,
                format=
                        "<level>{level}</level>: "  # 等级
                        "<cyan>{module}</cyan>.<cyan>{function}</cyan>"  # 模块名.方法名
                        ":<cyan>{line}</cyan> | "  # 行号
                        "<level>{message}</level>",  # 日志内容
                )

def set_batchsize(
    data_loader: _DataLoader, batch_size: int, encode=False
) -> _DataLoader:
    """Effectively copy the data loader, but with a different batch size.

    The `encode` option is used to copy the dataloader to use before encoding.
    This will not drop the last minibatch whose size may not be the requested
    batch size, and will also not shuffle the data.
    """
    return _DataLoader(
        dataset=data_loader.dataset,
        batch_size=batch_size,
        shuffle=not encode,
        drop_last=False,
        num_workers=1 if encode else data_loader.num_workers,
        pin_memory=data_loader.pin_memory,
    )


def make_dataloader(
    rpkm: _np.ndarray,
    tnf: _np.ndarray,
    dbf:_np.ndarray,
    lengths: _np.ndarray,
    batchsize: int = 256,
    destroy: bool = False,
    cuda: bool = True,
) -> _DataLoader:
    """Create a DataLoader from RPKM, TNF and lengths.

    The dataloader is an object feeding minibatches of contigs to the VAE.
    The data are normalized versions of the input datasets.

    Inputs:
        rpkm: RPKM matrix (N_contigs x N_samples)
        tnf: TNF matrix (N_contigs x N_TNF)
        lengths: Numpy array of sequence length (N_contigs)
        batchsize: Starting size of minibatches for dataloader
        destroy: Mutate rpkm and tnf array in-place instead of making a copy.
        cuda: Pagelock memory of dataloader (use when using GPU acceleration)

    Outputs:
        DataLoader: An object feeding data to the VAE
    """

    if not isinstance(rpkm, _np.ndarray) or not isinstance(tnf, _np.ndarray) or not isinstance(dbf,_np.ndarray):
        raise ValueError("TNF and RPKM must be Numpy arrays")

    if batchsize < 1:
        raise ValueError(f"Batch size must be minimum 1, not {batchsize}")

    if len(rpkm) != len(tnf) or len(tnf) != len(lengths) or len(dbf) != len(lengths):
        raise ValueError("Lengths of RPKM, TNF and lengths arrays must be the same")

    if not (rpkm.dtype == tnf.dtype== dbf.dtype == _np.float32):
        raise ValueError("TNF and RPKM must be Numpy arrays of dtype float32")

    # Copy if not destroy - this way we can have all following operations in-place
    # for simplicity
    if not destroy:
        rpkm = rpkm.copy()
        tnf = tnf.copy()
        dbf = dbf.copy()

    # Normalize samples to have same depth
    # 使得样本具有同样的深度，不会受到测序深度的影响
    sample_depths_sum = rpkm.sum(axis=0)
    if _np.any(sample_depths_sum == 0):
        raise ValueError(
            "One or more samples have zero depth in all sequences, so cannot be depth normalized"
        )
    rpkm *= 1_000_000 / sample_depths_sum
    total_abundance = rpkm.sum(axis=1)

    # Normalize rpkm to sum to 1

    n_samples = rpkm.shape[1]
    zero_total_abundance = total_abundance == 0
    rpkm[zero_total_abundance] = 1 / n_samples
    nonzero_total_abundance = total_abundance.copy()
    nonzero_total_abundance[zero_total_abundance] = 1.0
    rpkm /= nonzero_total_abundance.reshape((-1, 1))

    # Normalize TNF and total abundance to make SSE loss work better
    total_abundance = _np.log(total_abundance.clip(min=0.001))
    _vambtools.zscore(total_abundance, inplace=True)
    _vambtools.zscore(tnf, axis=0, inplace=True)
    _vambtools.zscore(dbf, axis=0, inplace=True) # dbf
    total_abundance.shape = (len(total_abundance), 1)

    # Create weights
    lengths = (lengths).astype(_np.float32)
    weights = _np.log(lengths).astype(_np.float32) - 5.0
    weights[weights < 2.0] = 2.0
    weights *= len(weights) / weights.sum()
    weights.shape = (len(weights), 1)

    ### Create final tensors and dataloader ###
    depthstensor = _torch.from_numpy(rpkm)  # this is a no-copy operation
    tnftensor = _torch.from_numpy(tnf)
    dbftensor = _torch.from_numpy(dbf)

    total_abundance_tensor = _torch.from_numpy(total_abundance)
    weightstensor = _torch.from_numpy(weights)
    n_workers = 4 if cuda else 1
    dataset = _TensorDataset(
        depthstensor, tnftensor, dbftensor, total_abundance_tensor, weightstensor
    )
    dataloader = _DataLoader(
        dataset=dataset,
        batch_size=batchsize,
        drop_last=False,
        shuffle=True,
        num_workers=n_workers,
        pin_memory=cuda,
    )

    return dataloader


class VAE(_nn.Module):
    """Variational autoencoder, subclass of torch.nn.Module.

    Instantiate with:
        nsamples: Number of samples in abundance matrix
        nhiddens: list of n_neurons in the hidden layers [None=Auto]
        nlatent: Number of neurons in the latent layer [32]
        alpha: Approximate starting TNF/(CE+TNF) ratio in loss. [None = Auto]
        beta: Multiply KLD by the inverse of this value [200]
        dropout: Probability of dropout on forward pass [0.2]
        cuda: Use CUDA (GPU accelerated training) [False]

    vae.trainmodel(dataloader, nepochs batchsteps, lrate, modelfile)
        Trains the model, returning None

    vae.encode(self, data_loader):
        Encodes the data in the data loader and returns the encoded matrix.

    If alpha or dropout is None and there is only one sample, they are set to
    0.99 and 0.0, respectively
    """

    def __init__(
        self,
        nsamples: int,
        nhiddens: Optional[list[int]] = None,
        nlatent: int = 32,
        alpha: Optional[float] = None,
        beta: float = 200.0,
        dropout: Optional[float] = 0.2,
        cuda: bool = True,
        seed: int = 0,
    ):
        if nlatent < 1:
            raise ValueError(f"Minimum 1 latent neuron, not {nlatent}")

        if nsamples < 1:
            raise ValueError(f"nsamples must be > 0, not {nsamples}")

        # If only 1 sample, we weigh alpha and nhiddens differently
        if alpha is None:
            alpha = 0.15 if nsamples > 1 else 0.50
        if nhiddens is None:
            nhiddens = [512, 512] if nsamples > 1 else [256, 256]

        if dropout is None:
            dropout = 0.2 if nsamples > 1 else 0.0

        if any(i < 1 for i in nhiddens):
            raise ValueError(f"Minimum 1 neuron per layer, not {min(nhiddens)}")

        if beta <= 0:
            raise ValueError(f"beta must be > 0, not {beta}")

        if not (0 < alpha < 1):
            raise ValueError(f"alpha must be 0 < alpha < 1, not {alpha}")

        if not (0 <= dropout < 1):
            raise ValueError(f"dropout must be 0 <= dropout < 1, not {dropout}")

        _torch.manual_seed(seed)
        self.rng = _torch.Generator()
        self.rng.manual_seed(seed)
        super(VAE, self).__init__()

        # Initialize simple attributes
        self.usecuda = cuda
        self.nsamples = nsamples
        self.ntnf = 103
        self.ndbf = 768
        self.alpha = alpha
        self.beta = beta
        self.nhiddens = nhiddens # [512, 512] if nsamples > 1 else [256, 256]
        self.nlatent = nlatent # 32
        self.dropout = dropout

        # Initialize lists for holding hidden layers
        self.encoderlayers = _nn.ModuleList()
        self.encodernorms = _nn.ModuleList()
        self.decoderlayers = _nn.ModuleList()
        self.decodernorms = _nn.ModuleList()

        # Add all other hidden layers
        for nin, nout in zip(
            # + 1 for the total abundance
            [self.nsamples + self.ntnf +self.ndbf + 1] + self.nhiddens,
            self.nhiddens,
        ):
            self.encoderlayers.append(_nn.Linear(nin, nout))
            self.encodernorms.append(_nn.BatchNorm1d(nout))

        # Latent layers
        self.mu = _nn.Linear(self.nhiddens[-1], self.nlatent)

        # Add first decoding layer
        for nin, nout in zip([self.nlatent] + self.nhiddens[::-1], self.nhiddens[::-1]):
            self.decoderlayers.append(_nn.Linear(nin, nout))
            self.decodernorms.append(_nn.BatchNorm1d(nout))

        # Reconstruction (output) layer. + 1 for the total abundance
        self.outputlayer = _nn.Linear(self.nhiddens[0], self.nsamples + self.ntnf + self.ndbf + 1)

        # Activation functions
        self.relu = _nn.LeakyReLU()
        self.softplus = _nn.Softplus()
        self.dropoutlayer = _nn.Dropout(p=self.dropout)

        if cuda:
            self.cuda()

    def _encode(self, tensor: Tensor) -> Tensor:
        tensors = list()

        # Hidden layers
        # 输入在经过每一层时，既被非线性转换（通过 Leaky ReLU），又被归一化处理（通过 BatchNorm），并且还应用了 Dropout 来减少过拟合风险。
        for encoderlayer, encodernorm in zip(self.encoderlayers, self.encodernorms):
            tensor = encodernorm(self.dropoutlayer(self.relu(encoderlayer(tensor))))
            tensors.append(tensor)

        # Latent layers
        mu = self.mu(tensor)

        # Note: We ought to also compute logsigma here, but we had a bug in the original
        # implementation of Vamb where logsigma was fixed to zero, so we just remove it.

        return mu

    # sample with gaussian noise
    def reparameterize(self, mu: Tensor) -> Tensor:
        epsilon = _torch.randn(mu.size(0), mu.size(1))

        if self.usecuda:
            epsilon = epsilon.cuda()

        epsilon.requires_grad = True

        latent = mu + epsilon

        return latent

    def _decode(self, tensor: Tensor) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        tensors = list()

        for decoderlayer, decodernorm in zip(self.decoderlayers, self.decodernorms):
            tensor = decodernorm(self.dropoutlayer(self.relu(decoderlayer(tensor))))
            tensors.append(tensor)

        reconstruction = self.outputlayer(tensor)

        # Decompose reconstruction to depths and tnf signal
        depths_out = reconstruction.narrow(1, 0, self.nsamples)
        tnf_out = reconstruction.narrow(1, self.nsamples, self.ntnf)
        dbf_out = reconstruction.narrow(1, self.nsamples, self.ndbf)
        abundance_out = reconstruction.narrow(1, self.nsamples + self.ntnf + self.ndbf, 1)

        depths_out = _softmax(depths_out, dim=1)

        return depths_out, tnf_out, dbf_out, abundance_out

    def forward(
        self, depths: Tensor, tnf: Tensor, dbf:Tensor, abundance: Tensor
    ) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
        tensor = _torch.cat((depths, tnf, dbf, abundance), 1)
        mu = self._encode(tensor)
        latent = self.reparameterize(mu)
        depths_out, tnf_out, dbf_out, abundance_out = self._decode(latent)

        return depths_out, tnf_out, dbf_out, abundance_out, mu

    def calc_loss(
        self,
        depths_in: Tensor,
        depths_out: Tensor,
        tnf_in: Tensor,
        tnf_out: Tensor,
        dbf_in:Tensor,
        dbf_out:Tensor,
        abundance_in: Tensor,
        abundance_out: Tensor,
        # 丰度描述的是样本中某种生物体或者基因的相对数量，
        # 而深度描述的是测序实验中数据的量或覆盖的程度
        mu: Tensor,
        weights: Tensor,
    ) -> tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
        ab_sse = (abundance_out - abundance_in).pow(2).sum(dim=1) 
        # Add 1e-9 to depths_out to avoid numerical instability.
        ce = -((depths_out + 1e-9).log() * depths_in).sum(dim=1) # 深度的交叉熵损失函数
        sse = (tnf_out - tnf_in).pow(2).sum(dim=1) # tnf的误差平方和
        kld = 0.5 * (mu.pow(2)).sum(dim=1) # KL散度
        dbe = (dbf_out-dbf_in).pow(2).sum(dim=1) # dbf也定义为误差平方和
        # Avoid having the denominator be zero
        # 单样本，ce_weight的权重设为0
        if self.nsamples == 1:
            ce_weight = 0.0
        else:
            ce_weight = ((1 - self.alpha) * (self.nsamples - 1)) / (
                self.nsamples * _log(self.nsamples)
            )

        ab_sse_weight = (1 - self.alpha) * (1 / self.nsamples)
        sse_weight = self.alpha / self.ntnf
        dbe_weight= self.alpha / self.ndbf
        kld_weight = 1 / (self.nlatent * self.beta)

        weighed_ab = ab_sse * ab_sse_weight
        weighed_ce = ce * ce_weight
        weighed_sse = sse * sse_weight
        weighted_dbe = dbe *dbe_weight
        weighed_kld = kld * kld_weight

        reconstruction_loss = weighed_ce + weighed_ab + weighed_sse + weighted_dbe
        loss = (reconstruction_loss + weighed_kld) * weights

        return (
            loss.mean(),
            weighed_ab.mean(),
            weighed_ce.mean(),
            weighed_sse.mean(),
            weighted_dbe.mean(),
            weighed_kld.mean(),
        )
###########################################################################
    def trainepoch(
        self,
        data_loader: _DataLoader,
        epoch: int,
        optimizer,
        batchsteps: list[int],
    ) -> _DataLoader[tuple[Tensor, Tensor, Tensor,Tensor]]:
        self.train()

        epoch_loss = 0.0
        epoch_kldloss = 0.0
        epoch_sseloss = 0.0
        epoch_dbeloss = 0.0
        epoch_celoss = 0.0
        epoch_absseloss = 0.0

        if epoch in batchsteps:
            data_loader = set_batchsize(
                data_loader, data_loader.batch_size * 2  # type:ignore
            )

        for depths_in, tnf_in, dbf_in, abundance_in, weights in data_loader:
            depths_in.requires_grad = True
            tnf_in.requires_grad = True
            dbf_in.requires_grad = True
            abundance_in.requires_grad = True

            if self.usecuda:
                depths_in = depths_in.cuda()
                tnf_in = tnf_in.cuda()
                dbf_in = dbf_in.cuda()
                abundance_in = abundance_in.cuda()
                weights = weights.cuda()

            optimizer.zero_grad()

            depths_out, tnf_out, dbf_out, abundance_out, mu = self(
                depths_in, tnf_in, dbf_in, abundance_in
            )

            loss, ab_sse, ce, sse, dbe, kld = self.calc_loss(
                depths_in,
                depths_out,
                tnf_in,
                tnf_out,
                dbf_in,
                dbf_out,
                abundance_in,
                abundance_out,
                mu,
                weights,
            )

            loss.backward()
            optimizer.step()

            epoch_loss += loss.data.item()
            epoch_kldloss += kld.data.item()
            epoch_sseloss += sse.data.item()
            epoch_dbeloss += dbe.data.item()
            epoch_celoss += ce.data.item()
            epoch_absseloss += ab_sse.data.item()
        
        logger.info(
            "\tEpoch: {:>3}  Loss: {:.5e}  CE: {:.5e}  AB: {:.5e}  SSE: {:.5e} DBE:{:.5e}  KLD: {:.5e}  Batchsize: {:>4}".format(
                epoch + 1,
                epoch_loss / len(data_loader),
                epoch_celoss / len(data_loader),
                epoch_absseloss / len(data_loader),
                epoch_sseloss / len(data_loader),
                epoch_dbeloss / len(data_loader),
                epoch_kldloss / len(data_loader),
                data_loader.batch_size,
            )
        )

        self.eval()
        return data_loader

    def encode(self, data_loader) -> _np.ndarray:
        """Encode a data loader to a latent representation with VAE

        Input: data_loader: As generated by train_vae

        Output: A (n_contigs x n_latent) Numpy array of latent repr.
        """
        logger.info("Encoding to latent representation.........")
        self.eval()

        new_data_loader = set_batchsize(
            data_loader, data_loader.batch_size, encode=True
        )

        depths_array, _, _, _, _ = data_loader.dataset.tensors
        length = len(depths_array)

        # We make a Numpy array instead of a Torch array because, if we create
        # a Torch array, then convert it to Numpy, Numpy will believe it doesn't
        # own the memory block, and array resizes will not be permitted.
        latent = _np.empty((length, self.nlatent), dtype=_np.float32)

        row = 0
        with _torch.no_grad():
            for depths, tnf, dbf, ab, _ in new_data_loader:
                # Move input to GPU if requested
                if self.usecuda:
                    depths = depths.cuda()
                    tnf = tnf.cuda()
                    dbf = dbf.cuda()
                    ab = ab.cuda()

                # Evaluate
                _, _, _, _, mu = self(depths, tnf, dbf, ab)

                if self.usecuda:
                    mu = mu.cpu()

                latent[row : row + len(mu)] = mu
                row += len(mu)

        assert row == length
        elapsed = round(time.time() - begintime, 2)
        logger.info(f"Trained VAE and encoded in {elapsed} seconds.\n")
        return latent

    def save(self, filehandle):
        """Saves the VAE to a path or binary opened file. Load with VAE.load

        Input: Path or binary opened filehandle
        Output: None
        """
        state = {
            "nsamples": self.nsamples,
            "alpha": self.alpha,
            "beta": self.beta,
            "dropout": self.dropout,
            "nhiddens": self.nhiddens,
            "nlatent": self.nlatent,
            "state": self.state_dict(),
        }

        _torch.save(state, filehandle)

    @classmethod
    def load(
        cls, path: Union[IO[bytes], str], cuda: bool = False, evaluate: bool = True
    ):
        """Instantiates a VAE from a model file.

        Inputs:
            path: Path to model file as created by functions VAE.save or
                  VAE.trainmodel.
            cuda: If network should work on GPU [False]
            evaluate: Return network in evaluation mode [True]

        Output: VAE with weights and parameters matching the saved network.
        """

        # Forcably load to CPU even if model was saves as GPU model
        dictionary = _torch.load(path, map_location=lambda storage, loc: storage)

        nsamples = dictionary["nsamples"]
        alpha = dictionary["alpha"]
        beta = dictionary["beta"]
        dropout = dictionary["dropout"]
        nhiddens = dictionary["nhiddens"]
        nlatent = dictionary["nlatent"]
        state = dictionary["state"]

        vae = cls(nsamples, nhiddens, nlatent, alpha, beta, dropout, cuda)
        vae.load_state_dict(state)

        if cuda:
            vae.cuda()

        if evaluate:
            vae.eval()

        return vae

    def trainmodel(
        self,
        dataloader: _DataLoader[tuple[Tensor, Tensor, Tensor, Tensor]], # 增加了一个Tensor
        nepochs: int = 200,
        lrate: float = 1e-3,
        batchsteps: Optional[list[int]] = [25, 75, 150],
        modelfile: Union[None, str, Path, IO[bytes]] = None,
    ):
        """Train the autoencoder from depths array and tnf array.

        Inputs:
            dataloader: DataLoader made by make_dataloader
            nepochs: Train for this many epochs before encoding [500]
            lrate: Starting learning rate for the optimizer [0.001]
            batchsteps: None or double batchsize at these epochs [25, 75, 150, 300]
            modelfile: Save models to this file if not None [None]

        Output: None
        """

        if lrate < 0:
            raise ValueError(f"Learning rate must be positive, not {lrate}")
        if nepochs < 1:
            raise ValueError("Minimum 1 epoch, not {nepochs}")

        if batchsteps is None:
            batchsteps_set: set[int] = set()
        else:
            # First collect to list in order to allow all element types, then check that
            # they are integers
            batchsteps = list(batchsteps)
            if not all(isinstance(i, int) for i in batchsteps):
                raise ValueError("All elements of batchsteps must be integers")
            if max(batchsteps, default=0) >= nepochs:
                raise ValueError("Max batchsteps must not equal or exceed nepochs")
            batchsteps_set = set(batchsteps)

        # Get number of features
        # Following line is un-inferrable due to typing problems with DataLoader
        ncontigs, nsamples = dataloader.dataset.tensors[0].shape  # type: ignore
        optimizer = _Adam(self.parameters(), lr=lrate)

        global begintime
        begintime = time.time()
        logger.info("-----------creating and training VAE-----------")
        logger.info("\tNetwork properties:")
        logger.info(f"\tCUDA: {self.usecuda}")
        logger.info(f"\tAlpha: {self.alpha}")
        logger.info(f"\tBeta: {self.beta}")
        logger.info(f"\tDropout: {self.dropout}")
        logger.info(f"\tN hidden: {', '.join(map(str, self.nhiddens))}")
        logger.info(f"\tN latent: {self.nlatent}")
        logger.info("\tTraining properties:")
        logger.info(f"\tN epochs: {nepochs}")
        logger.info(f"\tStarting batch size: {dataloader.batch_size}")
        batchsteps_string = (
            ", ".join(map(str, sorted(batchsteps_set))) if batchsteps_set else "None"
        )
        logger.info(f"\tBatchsteps: {batchsteps_string}")
        logger.info(f"\tLearning rate: {lrate}")
        logger.info(f"\tN sequences: {ncontigs}")
        logger.info(f"\tN samples: {nsamples}")
        logger.info("------------------trainepoch------------------")

        # Train
        for epoch in range(nepochs):
            dataloader = self.trainepoch(
                dataloader, epoch, optimizer, sorted(batchsteps_set)
            )

        # Save weights - Lord forgive me, for I have sinned when catching all exceptions
        if modelfile is not None:
            try:
                self.save(modelfile)
            except:
                pass

        return None
