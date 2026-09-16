#!/usr/bin/env python

"""
This module provides plotting functions for AnnData objects utilizing plotly.
"""

import itertools
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import scipy
from anndata import AnnData
from plotly.subplots import make_subplots

ERR_MSG_KEY_NOT_FOUND = "Could not find {} in {}."
ERR_MSG_UNEXPECTED_CATEGORIES = (
    "The following categories were not found in specified groups: {}\\n" "Specified groups: {}"
)


def qc_metrics(
    adata: AnnData,
    qc_vars: str | list | None = None,
    ncols: int = 3,
    *,
    width: int | None = None,
    height: int | None = None,
    bargap: int = 0,
    quantile: float = 1.0,
    template: str | None = None,
    wspace: float = 0.075,
    hspace: float = 0.1,
    return_fig: bool = False,
) -> go.Figure | None:
    """
    Calculate and create histograms of QC metrics:
        N counts per cell, N genes per cell, % of mitochondrial expression per cell or other qc_vars

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    qc_vars : str | list | None
        List of qc_vars specified, same as in sc.pp.calculate_metrics.
    ncols : int
        Maximum n of columns in subplots layout.
    width : int | None
        Width of the figure in pixels.
    height : int | None
        Height of the figure in pixels.
    bargap: int
        Gap between bars in histograms.
    quantile: float
        PLot data only of the first quantile (if specified).
    template : str | dict
        Plotly template name or template dict.
    wspace : float
        Spacing between subplot columns.
    hspace : float
        Spacing between subplot rows.
    return_fig : bool
        If True, return the figure instead of displaying it.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    if template is None:
        template = pio.templates.default

    if qc_vars:
        qc_vars = [qc_vars] if isinstance(qc_vars, str) else list(qc_vars)

    num_plots = 2 + len(qc_vars) if qc_vars else 2
    cols = min(num_plots, ncols)
    rows = int(np.ceil(num_plots / cols))

    subplot_titles = (
        ["N of UMI per cell", "N of genes per cell"] + [f"% of {qc_var} expression per cell" for qc_var in qc_vars]
        if qc_vars
        else ["N of UMI per cell", "N of genes per cell"]
    )

    fig = make_subplots(
        rows=rows, cols=cols, horizontal_spacing=wspace, vertical_spacing=hspace, subplot_titles=subplot_titles
    )

    all_metrics = (
        ["total_counts", "n_genes_by_counts"] + [f"pct_counts_{qc_var}" for qc_var in qc_vars]
        if qc_vars
        else ["total_counts", "n_genes_by_counts"]
    )

    for c, metric in enumerate(all_metrics):
        row = (c // cols) + 1
        col = (c % cols) + 1
        threshold = np.percentile(adata.obs[metric], quantile * 100)
        fig.add_trace(go.Histogram(x=adata.obs[metric][adata.obs[metric] <= threshold]), row=row, col=col)

    for row in range(1, rows + 1):
        for col in range(1, cols + 1):
            fig.update_yaxes(showticklabels=False, ticks="", showline=False, row=row, col=col)

    fig.update_layout(showlegend=False, template=template, width=width, height=height)
    fig.update_layout(bargap=bargap)
    fig.update_annotations(font={"family": "Serif"})

    if return_fig:
        return fig

    fig.show()
    return None


def _get_basis(adata: AnnData, basis: str) -> np.ndarray:
    """
    Get array for basis from anndata. Just tries to add 'X_'.
    """
    additional_message = (
        " Run sc.pp.pca first."
        if basis in ("pca", "X_pca")
        else " Run sc.tl.umap first." if basis in ("umap", "X_umap") else ""
    )
    if basis in adata.obsm:
        return adata.obsm[basis]
    if f"X_{basis}" in adata.obsm:
        return adata.obsm[f"X_{basis}"]
    raise KeyError(ERR_MSG_KEY_NOT_FOUND.format(f"{basis} or X_{basis}", ".obsm") + additional_message)


def _prepare_dimension_dataframes(adata, basis, dimensions, last_color_col, groups):
    """
    Prepare DataFrames for each dimension pair.
    """
    bas = basis[2:].upper() if basis.startswith("X_") else basis.upper()
    dfs_for_plot = []

    for dim_pair in dimensions:
        new_cols_df = pd.DataFrame(
            _get_basis(adata, basis)[:, dim_pair],
            columns=[f"{bas}{dim_pair[0]+1}", f"{bas}{dim_pair[1]+1}"],
        )
        df_for_plot = pd.concat(
            [adata.obs.reset_index(), new_cols_df],
            axis=1,
        ).set_index(adata.obs_names)
        for col in df_for_plot.select_dtypes(include="category").columns:
            df_for_plot[col] = df_for_plot[col].cat.add_categories("NA").fillna("NA")

        if groups:
            df_for_plot[last_color_col] = df_for_plot[last_color_col].where(
                df_for_plot[last_color_col].isin(groups),
                "NA",
            )

        dfs_for_plot.append(df_for_plot)

    return dfs_for_plot, bas


def _add_color_data(adata, dfs_for_plot, color):
    """
    Add color data to DataFrames.
    """
    for df_for_plot in dfs_for_plot:
        if color == [""]:
            df_for_plot[""] = pd.Categorical([""] * len(df_for_plot))
        else:
            for color_col in color:
                if color_col not in adata.obs.columns:
                    if color_col in adata.var_names:
                        df_for_plot[color_col] = adata[:, adata.var_names == color_col].X.toarray().T[0]
                    else:
                        raise KeyError(ERR_MSG_KEY_NOT_FOUND.format(color_col, ".var_names or .obs.columns"))

    return color


def _add_trace_to_figure(
    fig, trace, row, col, counter, color_col, dim_pair, *, group_legends: bool, shared_coloraxes: bool
):
    """
    Add a trace to the figure with proper formatting.
    """
    if not group_legends:
        trace.update(legend=f"legend{counter + 1}")

    fig.add_trace(trace, row=row, col=col)

    if group_legends:
        legendgroup_name = f"{color_col} {dim_pair}" if dim_pair != (0, 1) else color_col
        fig.update_traces(
            legendgrouptitle=dict(text=legendgroup_name),
            legendgroup=legendgroup_name,
            row=row,
            col=col,
        )
    if not shared_coloraxes:
        fig.update_traces(marker=dict(coloraxis=f"coloraxis{counter+1}"))
        fig.update_layout({f"coloraxis{counter+1}": {"showscale": False}})


def _calculate_centroids(df_for_plot, color_col, x_col, y_col):
    """
    Calculate centroids of each group for a given color_col.
    """
    return df_for_plot.groupby(color_col, observed=False)[[x_col, y_col]].mean()


def _add_annotations(fig, centroids, row, col):
    """
    Add annotations (category labels) at centroids.
    """
    for cat, coords in centroids.iterrows():
        if cat != "NA":
            fig.add_annotation(
                x=coords.iloc[0],
                y=coords.iloc[1],
                text=str(cat),
                showarrow=False,
                xref=f"x{col}",
                yref=f"y{row}",
                name="annotations",
            )


def _update_annotations(fig, annotations_font, annotations_outline_width, subtitles_font):
    annotations_font = annotations_font or {"family": "Serif"}
    subtitles_font = subtitles_font or {"size": 16, "family": "Serif"}
    if annotations_outline_width:
        annotations_font["shadow"] = (
            f"{annotations_outline_width}px {annotations_outline_width}px 0 white, "
            f"-{annotations_outline_width}px {annotations_outline_width}px 0 white, "
            f"{annotations_outline_width}px -{annotations_outline_width}px 0 white, "
            f"-{annotations_outline_width}px -{annotations_outline_width}px 0 white, "
            f"0 {annotations_outline_width}px 0 white, "
            f"0 -{annotations_outline_width}px 0 white, "
            f"{annotations_outline_width}px 0 0 white, "
            f"-{annotations_outline_width}px 0 0 white"
        )
    fig.update_annotations(font=subtitles_font)
    fig.update_annotations(selector={"name": "annotations"}, font=annotations_font)


def _update_legends(fig, rows, cols, num_plots, hspace, wspace, showcoloraxes, shared_coloraxes, shared_legends, cmap):
    if shared_coloraxes:
        fig.layout.coloraxis.colorbar.x = 1.05
        fig.layout.legend.x = 1.1
        fig.layout.coloraxis.colorbar.xref = "paper"
        fig.layout.coloraxis.colorbar.yref = "paper"
        fig.layout.coloraxis.colorbar.y = 0.5
        fig.layout.coloraxis.colorbar.thickness = 10
        fig.layout.coloraxis.colorbar.len = 1 / rows
        if cmap:
            fig.layout.coloraxis.colorscale = cmap
    else:
        for i in range(1, num_plots + 1):
            coloraxis = f"coloraxis{i}"

            # Determine row and column for the current subplot
            row = (i - 1) // cols  # 0-indexed row position
            col = (i - 1) % cols  # 0-indexed column position

            # Calculate x position based on the column, spacing, and total number of columns
            x_position = (col + 1) / cols - (wspace * (cols - col - 1) / cols)

            # Calculate y position based on the row, spacing, and total number of rows
            y_position = 1 - (row + 0.5) / rows * (1 + hspace) + (hspace / 2)

            # Set the color axis properties
            fig.layout[coloraxis].showscale = showcoloraxes
            fig.layout[coloraxis].colorbar.x = x_position
            fig.layout[coloraxis].colorbar.y = y_position

            # Set other properties, such as thickness and length
            fig.layout[coloraxis].colorbar.thickness = 10
            fig.layout[coloraxis].colorbar.len = 0.75 / rows
            if i <= len(fig.data):
                fig.data[i - 1].update(marker=dict(coloraxis=coloraxis))
            if cmap:
                fig.layout[f"coloraxis{i}"].colorscale = cmap

    if shared_legends:
        fig.update_layout(legend=dict(x=1.05, groupclick="toggleitem"))
        if shared_coloraxes:
            fig.update_layout(legend=dict(x=1.12, groupclick="toggleitem"))

    else:
        for i in range(1, num_plots + 1):
            legend = f"legend{i}"
            row = (i - 1) // cols  # 0-indexed row position
            col = (i - 1) % cols  # 0-indexed column position
            x_position = (col + 1) / cols - (wspace * (cols - col - 1) / cols)
            y_position = 1 - (row + 0.5) / rows * (1 + hspace) + (hspace / 2)
            fig.update_layout(
                {
                    legend: dict(
                        x=x_position,
                        y=y_position,
                        yanchor="middle",
                    )
                }
            )


def _modify_cmap(cmap, zero_color):
    if zero_color:
        if cmap:
            prev_cmap = pc.get_colorscale(cmap) if isinstance(cmap, str) else cmap
        else:
            prev_cmap = pio.templates[pio.templates.default].layout.colorscale.sequential
        return [*[[0, zero_color], [0.00001, prev_cmap[0][1]]], *prev_cmap[1:]]
    return cmap


def embedding(
    adata: AnnData,
    basis: str,
    *,
    marker_size: float | None = 3,
    marker_edgewidth: float | None = None,
    marker_edgecolor: str | None = None,
    template: str | None = None,
    dimensions: tuple[int, int] | list[tuple[int, int]] = (0, 1),
    color: str | list[str] | None = None,
    groups: list[str] | None = None,
    annotations: bool = False,
    annotations_font: dict | None = None,
    annotations_outline_width: int | None = None,
    ncols: int = 3,
    wspace: float = 0.1,
    hspace: float = 0.15,
    opacity: float = 1,
    hover_name: str | pd.Series | None = None,
    hover_data: str | list[str] | pd.Series | dict | None = None,
    shared_axes: bool | str = "all",
    shared_coloraxes: bool = False,
    shared_legends: bool = False,
    width: int | None = None,
    height: int | None = None,
    subtitles: str | list[str] | None = None,
    subtitles_font: dict | None = None,
    cmap: str | list[list[int]] | None = None,
    zero_color: str | None = None,
    na_color: str = "lightgray",
    showlegend: bool | None = None,
    showcoloraxes: bool = True,
    return_fig: bool = False,
    _pca_annotate_variances: bool = False,
    **kwargs,
) -> go.Figure | None:
    """
    Create embedding plots with multiple dimensions and colors.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    basis : str
        Basis for dimensional reduction (e.g., 'X_umap', 'X_tsne', 'X_pca').
    marker_size : float
        Size of markers in scatter plots.
    marker_edgewidth : float | None
        Edge width of markers in scatter plots.
    marker_edgecolor : str | None
        Edge color of markers in scatter plots.
    template : str | dict
        Plotly template name or template dict.
    dimensions : tuple[int, int] | list[tuple[int, int]]
        Dimensions to plot, either single tuple or list of tuples.
    color : str | list[str] | None
        Column name(s) to use for coloring.
    groups: list[str] | None
        Show only these groups from the last color in scatter plots.
    annotations: bool
        Show each group name in their centroids in scatter plots.
    annotations_font: dict | None
        Font parameters of annotations.
    annotations_outline_width: int | None
        White outline for text (actually created using css shadows).
    ncols : int
        Maximum number of columns in subplot grid.
    wspace : float
        Spacing between subplot columns.
    hspace : float
        Spacing between subplot rows.
    opacity : float
        Opacity of markers (0 to 1).
    hover_name : str | int | pd.Series | np.array| None
        Hover_name for px.express, shows the data bold in the hover tooltip.
    hover_data : str | list[str | int] | pd.Series | np.array | dict | None
        Hover_data for px.express, shows the data in the hover tooltip.
    shared_axes: bool | str
        True or "columns": Share axes among subplots in the same column.
        "rows": Share axes among subplots in the same row
        "all": Share axes across all subplots in the grid.
        False: don't share axes.
    shared_coloraxes: bool
        Use the same coloraxes for all subplots.
    shared_legends: bool
        Display all legends on the right (using legend groups).
    width : int | None
        Width of figure in pixels.
    height : int | None
        Height of figure in pixels.
    subtitles : str | list[str] | None
        List of subtitles to use, color by default.
    subtitles_font : dict | None
        Font parameters of subtitles.
    cmap : str | None
        Colormap for continous variables.
    zero_color : str | None
        Color to use for 0 values while plotting continous variables.
    na_color : str
        Color to use for NA values.
    showlegend : bool
        Show all legends.
    showcoloraxes : bool
        Show all coloraxes.
    return_fig : bool
        If True, return the figure instead of displaying it.
    **kwargs
        Additional keyword arguments passed to fig.update_layout.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    template = template or pio.templates.default
    cmap = _modify_cmap(cmap, zero_color)

    # Prepare data
    if isinstance(dimensions, tuple):
        dimensions = [dimensions]
    if color is None:
        color = [""]
    else:
        color = [color] if isinstance(color, str) else list(color)

    last_color_col = color[-1]
    dfs_for_plot, bas = _prepare_dimension_dataframes(adata, basis, dimensions, last_color_col, groups)
    color = _add_color_data(adata, dfs_for_plot, color)

    # Calculate layout
    num_plots = len(color) * len(dimensions)
    cols = min(num_plots, ncols)
    rows = int(np.ceil(num_plots / cols))

    # Create figure
    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=(
            [subtitles]
            if isinstance(subtitles, str)
            else list(subtitles) if subtitles is not None else np.repeat(color, len(dimensions))
        ),
        horizontal_spacing=wspace,
        vertical_spacing=hspace,
        shared_xaxes=shared_axes,
        shared_yaxes=shared_axes,
    )

    # Add subplots
    counter = 0
    for color_col in color:
        for dim_pair, df_for_plot in zip(dimensions, dfs_for_plot):
            x_col = f"{bas}{dim_pair[0]+1}"
            y_col = f"{bas}{dim_pair[1]+1}"
            if (color_col in adata.obs) and (adata.obs[color_col].dtype.name == "category"):
                category_orders = {color_col: adata.obs[color_col].cat.categories}
            else:
                category_orders = None

            px_fig = px.scatter(
                df_for_plot,
                x=x_col,
                y=y_col,
                template=template,
                color=color_col,
                opacity=opacity,
                hover_name=df_for_plot.index if hover_name is None else hover_name,
                hover_data=hover_data,
                color_discrete_map={"NA": na_color},
                category_orders=category_orders,
            )

            for trace in px_fig["data"]:
                row = (counter // cols) + 1
                col = (counter % cols) + 1
                _add_trace_to_figure(
                    fig,
                    trace,
                    row,
                    col,
                    counter,
                    color_col,
                    dim_pair,
                    group_legends=shared_legends,
                    shared_coloraxes=shared_coloraxes,
                )

            if annotations and df_for_plot[color_col].dtype.name == "category":
                centroids = _calculate_centroids(df_for_plot, color_col, x_col, y_col)
                _add_annotations(fig, centroids, row, col)

            if _pca_annotate_variances:
                var_x = round(adata.uns["pca"]["variance_ratio"][dim_pair[0]] * 100, 2)
                var_y = round(adata.uns["pca"]["variance_ratio"][dim_pair[1]] * 100, 2)
                fig.update_xaxes(title_text=f"{x_col} ({var_x}%)", row=row, col=col)
                fig.update_yaxes(title_text=f"{y_col} ({var_y}%)", row=row, col=col)
            else:
                fig.update_xaxes(title_text=x_col, row=row, col=col)
                fig.update_yaxes(title_text=y_col, row=row, col=col)
            counter += 1

    # fig.layout.legend.x = 1.05

    _update_annotations(fig, annotations_font, annotations_outline_width, subtitles_font)

    _update_legends(fig, rows, cols, num_plots, hspace, wspace, showcoloraxes, shared_coloraxes, shared_legends, cmap)
    fig.update_yaxes(showticklabels=False, zeroline=False, ticks="")
    fig.update_xaxes(showticklabels=False, zeroline=False, ticks="")
    fig.update_traces(marker_size=marker_size, marker_line=dict(width=marker_edgewidth, color=marker_edgecolor))

    showlegend = (
        False if showlegend is None and annotations else True if showlegend is None and not annotations else showlegend
    )
    fig.update_layout(
        showlegend=showlegend, template=template, margin={"pad": 20}, width=width, height=height, **kwargs
    )

    if return_fig:
        return fig
    fig.show()
    return None


def pca(
    adata: AnnData,
    *,
    annotate_var_explained: bool = True,
    **kwargs,
) -> go.Figure | None:
    """
    Scatter plot in PCA coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    annotate_var_explained : bool
        Whether to annotate the explained variance on axis names
    **kwargs
        Additional keyword arguments passed to embedding.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    return embedding(adata, basis="pca", _pca_annotate_variances=annotate_var_explained, **kwargs)


def umap(adata: AnnData, **kwargs):
    """
    Scatter plot in UMAP basis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    **kwargs
        Additional keyword arguments passed to embedding.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    return embedding(adata, basis="umap", **kwargs)


def highly_variable_genes(
    adata: AnnData,
    *,
    log: bool = True,
    shared_axes: bool = False,
    opacity: float = 1,
    marker_size: float | None = 3,
    marker_edgewidth: float | None = None,
    marker_edgecolor: str | None = None,
    width: int | None = None,
    height: int | None = None,
    template: str | None = None,
    return_fig: bool = False,
):
    """
    Create scatter plots comparing normalized and non-normalized variances of genes
    against their mean expression, highlighting highly variable genes.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain 'means', 'variances', 'variances_norm',
        and 'highly_variable' columns in adata.var.
    log : bool
        If True, use logarithmic scale for both axes.
    shared_axes : bool
        If True, share the same scale across both plots.
    opacity : float
        Opacity of markers in scatter plots (0 to 1).
    marker_size : float
        Size of markers in scatter plots.
    marker_edgewidth : float | None
        Edge width of markers in scatter plots.
    marker_edgecolor : str | None
        Edge color of markers in scatter plots.
    width : int | None
        Width of the figure in pixels.
    height : int | None
        Height of the figure in pixels.
    template : str | dict
        Plotly template name or template dict.
    return_fig : bool
        If True, return the figure instead of displaying it.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    template = template or pio.templates.default

    shared_axes = shared_axes if not shared_axes else "rows"
    fig = make_subplots(rows=1, cols=2, shared_xaxes=shared_axes, shared_yaxes=shared_axes)
    y = "variances_norm" if "variances_norm" in adata.var else "dispersions_norm"
    norm_var_px = px.scatter(
        adata.var,
        x="means",
        y=y,
        color="highly_variable",
        hover_name=adata.var_names,
        opacity=opacity,
        category_orders={"highly_variable": [True, False]},
        color_discrete_sequence=[
            pio.templates[pio.templates.default].layout.colorway[1],
            pio.templates[pio.templates.default].layout.colorway[0],
        ],
    )
    y = "variances" if "variances" in adata.var else "dispersions"
    var_px = px.scatter(
        adata.var,
        x="means",
        y=y,
        color="highly_variable",
        hover_name=adata.var_names,
        opacity=opacity,
        category_orders={"highly_variable": [True, False]},
        color_discrete_sequence=[
            pio.templates[pio.templates.default].layout.colorway[1],
            pio.templates[pio.templates.default].layout.colorway[0],
        ],
    )
    for trace in norm_var_px["data"]:
        fig.add_trace(trace, row=1, col=1)
        fig.update_traces(showlegend=False)
    for trace in var_px["data"]:
        fig.add_trace(trace, row=1, col=2)

    fig.for_each_trace(lambda t: t.update(name={"False": "other", "True": "highly variable"}[t.name]))
    fig.update_layout(
        width=width, height=height, legend_title_text="", template=template, title="Highly variable genes"
    )
    fig.update_traces(marker_size=marker_size, marker_line=dict(width=marker_edgewidth, color=marker_edgecolor))
    if log:
        fig.update_xaxes(type="log")
        fig.update_yaxes(type="log")
    fig.update_yaxes(title="variance (normalized)", row=1, col=1)
    fig.update_yaxes(title="variance (not normalized)", row=1, col=2)
    fig.update_xaxes(title="mean expression")
    if return_fig:
        return fig
    fig.show()
    return None


def pca_variance_ratio(
    adata: AnnData,
    n_pcs: int | None = None,
    *,
    log: bool = False,
    template: str | None = None,
    annotations: bool = True,
    return_fig: bool = False,
    **kwargs,
):
    """
    Create a scree plot showing explained variance ratio for principal components.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain PCA results in adata.uns['pca']
        with 'variance_ratio' key.
    n_pcs : int | None
        Number of principal components to plot. If None, all available
        components will be shown.
    log : bool
        Logarithmic y axes.
    template : str | dict
        Plotly template name or template dict.
    annotations : bool
        Create text annotations for each PC (True by default).
    return_fig : bool
        If True, return the figure instead of displaying it.
    **kwargs
        Additional keyword arguments passed to fig.update_layout.

    Returns
    -------
    go.Figure | None
        The figure object if return_fig is True, None otherwise.
    """
    template = template or pio.templates.default

    y = adata.uns["pca"]["variance_ratio"][:n_pcs]
    x = np.arange(1, len(y) + 1)
    plot_df = pd.DataFrame({"ranking": x, "explained variance": y, "PC": [f"PC{i}" for i in x]}).set_index("ranking")
    fig = px.line(plot_df, y="explained variance", markers=True, hover_name="PC")
    if log:
        fig.update_yaxes(type="log", dtick=1)
    if annotations:
        for c, val in enumerate(x):
            fig.add_annotation(
                x=val,
                y=np.log10(y[c]) if log else y[c],
                text=f"PC{val}",
                showarrow=False,
                xanchor="left",
                yanchor="bottom",
                textangle=90,
            )
    fig.update_xaxes(range=(0, len(y) + 1))
    fig.update_layout(yaxis_tickformat=",.0%", title=" PCA scree plot", template=template, **kwargs)
    if return_fig:
        return fig
    fig.show()
    return None


def _create_group_combinations(adata, groupby):
    """
    create group combinations if multiple groupby columns
    """
    group_combinations = []
    if len(groupby) > 1:
        unique_values = [adata.obs[col].unique() for col in groupby]
        combinations = list(itertools.product(*unique_values))
        group_combinations = ["_".join(f"{val}" for val in comb) for comb in combinations]
    else:
        group_combinations = adata.obs[groupby[0]].unique()
    return group_combinations


def _add_category_labels_annotations(fig, var_names, plot_df, dendrogram):
    if isinstance(var_names, dict):
        unique_cats = []
        current_pos = -0.5
        cat_positions = []
        cat_ranges = []
        for cat, genes in var_names.items():
            unique_cats.append(cat)
            n_genes = len(plot_df[plot_df["gene"].isin(genes)]["gene"].unique().tolist())
            start_pos = current_pos + 0.25
            end_pos = start_pos + n_genes - 0.25
            cat_position = (start_pos + end_pos) / 2

            cat_positions.append(cat_position)
            cat_ranges.append((start_pos, end_pos))
            current_pos = end_pos

        all_shapes = []

        yref = "y2" if dendrogram else "y"

        for cat, pos, (start_pos, end_pos) in zip(unique_cats, cat_positions, cat_ranges):
            # Add bracket lines
            all_shapes.extend(
                [
                    dict(
                        type="line",
                        xref="paper",
                        yref=yref,
                        x0=1.025,
                        x1=1.025,
                        y0=start_pos,
                        y1=end_pos,
                        line=dict(color="black", width=1),
                    ),
                    dict(
                        type="line",
                        xref="paper",
                        yref=yref,
                        x0=1,
                        x1=1.025,
                        y0=start_pos,
                        y1=start_pos,
                        line=dict(color="black", width=1),
                    ),
                    dict(
                        type="line",
                        xref="paper",
                        yref=yref,
                        x0=1,
                        x1=1.025,
                        y0=end_pos,
                        y1=end_pos,
                        line=dict(color="black", width=1),
                    ),
                ]
            )

            # Add rotated text label
            fig.add_annotation(
                x=1.05,
                y=pos,
                text=cat,
                textangle=90,
                xref="paper",
                yref=yref,
                showarrow=False,
                font=dict(size=12),
                xanchor="center",
                yanchor="middle",
            )
        fig.update_layout(shapes=all_shapes)


def _create_plot_data(adata, var_names, groupby):
    """
    Create plotting df for dotplot function.
    """
    group_combinations = _create_group_combinations(adata, groupby)
    genes_df = adata.var[["means"]].copy()
    n_cells_by_counts = (adata.X > 0).sum(axis=0)
    genes_df["n_cells_by_counts"] = (
        n_cells_by_counts.A1 if isinstance(n_cells_by_counts, np.matrix) else n_cells_by_counts
    )
    genes_df["pct_cells"] = genes_df["n_cells_by_counts"] / adata.n_obs * 100

    if isinstance(var_names, dict):
        genes_flat = []
        categories = []
        for category, genes in var_names.items():
            genes_flat.extend(genes)
            categories.extend([category] * len(genes))
    else:
        genes_flat = var_names
        categories = [None] * len(genes_flat)

    plot_data = []
    for group_combo in group_combinations:
        if len(groupby) > 1:
            # Create mask for multiple groupby columns
            mask = np.ones(len(adata), dtype=bool)
            for col, val in zip(groupby, group_combo.split("_"), strict=True):
                mask &= adata.obs[col] == val
        else:
            # Single groupby column
            mask = adata.obs[groupby[0]] == group_combo

        # Get subset for this group
        adata_group = adata[mask]

        # Calculate statistics for each gene in this group
        group_stats = pd.DataFrame(index=genes_flat)
        group_stats["means"] = (
            adata_group[:, genes_flat].X.mean(axis=0).A1
            if scipy.sparse.issparse(adata_group.X)
            else adata_group[:, genes_flat].X.mean(axis=0)
        )
        group_stats["n_cells_by_counts"] = (
            (adata_group[:, genes_flat].X > 0).sum(axis=0).A1
            if scipy.sparse.issparse(adata_group.X)
            else (adata_group[:, genes_flat].X > 0).sum(axis=0)
        )
        group_stats["pct_cells"] = group_stats["n_cells_by_counts"] / len(adata_group) * 100

        for gene in genes_flat:
            plot_data.append(
                {  # noqa: PERF401
                    "gene": gene,
                    "category": categories[genes_flat.index(gene)],
                    "group": str(group_combo),
                    "cells_fraction": group_stats.loc[gene, "pct_cells"],
                    "mean_expression": group_stats.loc[gene, "means"],
                }
            )

    return pd.DataFrame(plot_data)


def dotplot(
    adata: AnnData,
    var_names: list[str] | dict[str, str],
    groupby: str | list[str],
    *,
    dendrogram: bool = False,
    categories_order: list | None = None,
    size_max: int = 15,
    marker_edgewidth: float | None = 0.5,
    marker_edgecolor: str | None = "darkslategray",
    title: str = "",
    colorbar_title: str = "Mean expression",
    cmap: str | None = None,
    height: int | None = None,
    width: int | None = None,
    return_fig: bool = False,
    template: str | None = None,
    **kwargs,
):
    """
    Create a dot plot visualization with optional dendrogram support
    (use sc.tl.dendrogram beforehand to be able to plot dendrogram).

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    var_names : list | dict
        Variables (genes) to plot, should be in adata.var_names
    groupby : str | list
        Column name(s) to group by
    dendrogram : bool
        Whether to show dendrogram
    categories_order : list | None
        Custom category ordering
    size_max : int
        Maximum dot size
    marker_edgewidth : float | None
        Edge width of dots
    marker_edgecolor : str | None
        Edge color of dots
    title : str
        Title for the whole plot
    colorbar_title : str
        Colorbar title, "Mean expression" be default
    cmap : str | None
        Coloraxes colorscale (colormap) name
    height : int | None
        Plot height
    width : int | None
        Plot width
    return_fig : bool
        Whether to return the figure object
    template : str | None
        Plotly template name
    **kwargs : dict
        Additional arguments for fig.update_layout()

    Returns:
    --------
    plotly.graph_objects.Figure or None
        Figure object if return_fig is True, None otherwise
    """

    template = template or pio.templates.default

    groupby = [groupby] if isinstance(groupby, str) else list(groupby)
    if dendrogram:
        if "_".join(["dendrogram", *groupby]) not in adata.uns:
            raise KeyError(
                ERR_MSG_KEY_NOT_FOUND.format("_".join(["dendrogram", *groupby]), ".uns")
                + " Run sc.tl.dendrogram first."
            )
        dendrogram_data = adata.uns["_".join(["dendrogram", *groupby])]
    categories_order = dendrogram_data["categories_ordered"] if dendrogram else categories_order

    plot_df = _create_plot_data(adata, var_names, groupby)

    if categories_order:
        unmatched_categories = set(categories_order) - set(plot_df["group"].unique())
        if unmatched_categories:
            raise KeyError(ERR_MSG_UNEXPECTED_CATEGORIES.format(unmatched_categories, plot_df["group"].unique()))

    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=plot_df["gene"].unique(), ordered=True)
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=categories_order, ordered=True)
    plot_df = plot_df.sort_values(["group", "gene"])

    # Create plot
    scatter_fig = px.scatter(
        plot_df,
        x="group",
        y="gene",
        size="cells_fraction",
        color="mean_expression",
        custom_data=["cells_fraction", "mean_expression"],
        size_max=size_max,
        range_x=[-0.5, len(plot_df["group"].unique()) - 0.5],
        range_y=[-0.5, len(plot_df["gene"].unique()) - 0.5],
        width=width,
        height=height,
    )

    if not dendrogram:
        scatter_fig.update_layout(
            xaxis_title="",
            yaxis_title="",
            annotations=[
                dict(
                    text="Dot size — % cells",
                    xref="paper",
                    yref="paper",
                    x=1.075,
                    y=0.5,
                    xanchor="left",
                    yanchor="middle",
                    showarrow=False,
                    textangle=-90,  # Vertical text
                )
            ],
        )
        fig = scatter_fig
        fig.update_xaxes(type="category")

    else:
        dendrogram_height_ratio = 100 / scatter_fig.layout.height if scatter_fig.layout.height else 0.15
        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            row_heights=[dendrogram_height_ratio, 1 - dendrogram_height_ratio],
            vertical_spacing=0,
        )

        # Add dendrogram trace
        dendr_xs = dendrogram_data["dendrogram_info"]["icoord"]
        dendr_ys = dendrogram_data["dendrogram_info"]["dcoord"]
        for xs, ys in zip(dendr_xs, dendr_ys):
            fig.add_trace(
                go.Scatter(
                    x=xs, y=ys, mode="lines", hoverinfo="skip", line=dict(color="black", width=1), showlegend=False
                ),
                row=1,
                col=1,
            )

        for trace in scatter_fig.data:
            fig.add_trace(trace, row=2, col=1)

        factor = (np.max(dendr_xs) - np.min(dendr_xs)) / len(plot_df["group"].unique())

        fig.update_xaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[-0.5 * factor + np.min(dendr_xs), 0.5 * factor + np.max(dendr_xs)],
            showline=False,
            row=1,
            col=1,
        )
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, showline=False, row=1, col=1)
        fig.update_xaxes(title="", type="category", range=[-0.5, len(plot_df["group"].unique()) - 0.5], row=2, col=1)
        fig.update_yaxes(title="", range=[-0.5, len(plot_df["gene"].unique()) - 0.5], row=2, col=1)
        fig.update_layout(
            annotations=[
                dict(
                    text="Dot size — % cells",
                    xref="paper",
                    yref="paper",
                    x=1.075,
                    y=0.5,
                    xanchor="left",
                    yanchor="middle",
                    showarrow=False,
                    textangle=-90,  # Vertical text
                )
            ]
        )

    _add_category_labels_annotations(fig, var_names, plot_df, dendrogram)

    fig.update_coloraxes(
        reversescale=True,
        colorscale=cmap,
        colorbar_title=colorbar_title,
        colorbar_lenmode="pixels",
        colorbar_len=200,
        colorbar_thickness=20,
        colorbar_x=1.1,
    )

    fig.update_layout(height=height, width=width, title=title, template=template, **kwargs)

    # Customize hover template
    fig.update_traces(
        hovertemplate="<br>".join(
            [
                "Group: %{x}",
                "Gene: %{y}",
                "Cells expressing: %{customdata[0]:.1f}%",
                "Mean expression: %{customdata[1]:.2f}",
                "<extra></extra>",
            ]
        ),
        marker_line=dict(width=marker_edgewidth, color=marker_edgecolor),
    )

    if return_fig:
        return fig
    fig.show()
    return None


def _calculate_n_degs(deg_df: pd.DataFrame, logfc: float = 0.3, pval: float = 0.05) -> int:
    """
    deg_df:        table of DEGs, should include columns 'LFC' and 'pval_adj'
    logfc:     threshold for DEGs by absolute value logfc
    pval:      threshold for DEGs by adjusted p-value

    returns n of DEGs
    """
    return len(deg_df[(abs(deg_df["LFC"]) > logfc) & (deg_df["pval_adj"] < pval)])


def _process_df_for_volcano(
    deg_df: pd.DataFrame, logfc: float = 0.3, pval: float = 0.05, *, remove_outliers: bool = True
) -> pd.DataFrame:
    """
    rename columns, add columns '-log10 P-value' and 'DEG type',
    remove outliers using interquantile range method

    deg_df:            table of DEGs, should include columns 'LFC' and 'pval_adj'
    logfc:             threshold for DEGs by absolute value logfc
    pval:              threshold for DEGs by adjusted p-value
    remove_outliers:   flag to remove genes with logfc value >= quartile + 15 * interquartile range
    """
    deg_df = deg_df.copy()
    if remove_outliers:
        q1 = deg_df["LFC"].quantile(0.25)
        q3 = deg_df["LFC"].quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 15 * iqr
        upper_bound = q3 + 15 * iqr
        were_genes = len(deg_df)
        deg_df = deg_df[(deg_df["LFC"] >= lower_bound) & (deg_df["LFC"] <= upper_bound)]
        n_removed_genes = were_genes - len(deg_df)
        if n_removed_genes != 0:
            logging.warning("%s genes were removed as outliers", n_removed_genes)

    deg_df["DEG type"] = deg_df.apply(
        lambda x: (
            "Down-regulated genes"
            if ((x["LFC"] < -logfc) and (x["pval_adj"] < pval))
            else "Up-regulated genes" if ((x["LFC"] > logfc) and (x["pval_adj"] < pval)) else "Insignificant genes"
        ),
        axis=1,
    )
    deg_df["-log10 P-value"] = -np.log10(deg_df["pval_adj"])
    return deg_df.rename(columns={"LFC": "log Fold Change"})


def _plot_degs(
    deg_df: pd.DataFrame,
    logfc: float = 0.3,
    pval: float = 0.05,
    title: str | None = None,
    color_discrete_map=None,
    **kwargs,
):
    """
    volcano plot
    """
    color_discrete_map = color_discrete_map or {
        "Up-regulated genes": "seagreen",
        "Down-regulated genes": "pink",
        "Insignificant genes": "darkgray",
    }
    fig = px.scatter(
        deg_df,
        x="log Fold Change",
        y="-log10 P-value",
        color="DEG type",
        hover_name="Gene",
        hover_data=["pval_adj"],
        title=title,
        color_discrete_map=color_discrete_map,
        **kwargs,
    )
    fig.update_yaxes(range=[0, deg_df["-log10 P-value"].max() + 1])

    fig.add_shape(
        type="line",
        x0=logfc,
        x1=logfc,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color="lightgray", width=0.5),
    )
    fig.add_shape(
        type="line",
        x0=-logfc,
        x1=-logfc,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color="lightgray", width=0.5),
    )
    fig.add_shape(
        type="line",
        x0=0,
        x1=1,
        y0=-np.log10(pval),
        y1=-np.log10(pval),
        xref="paper",
        yref="y",
        line=dict(color="lightgray", width=0.5),
    )

    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
    )
    return fig


def volcano(
    degs_df: pd.DataFrame,
    logfc: float = 0.3,
    pval: float = 0.05,
    *,
    logfc_col: str = "logfoldchanges",
    padj_col: str = "pvals_adj",
    gene_col: str = "names",
    title: str | None = None,
    remove_outliers: bool = True,
    color_discrete_map: dict | None = None,
    opacity: float = 0.5,
    return_fig: bool = False,
    template: str | None = None,
    **kwargs,
):
    """
    Create a volcano visualization of differentially expressed genes

    Parameters:
    -----------
    degs_df : DataFrame
        Dataframe wirh DEGs (should have logfc, pvalue and gene name columns).
    logfc : float
        Log fold change threshold for significant genes.
    pval : float
        P-value threshold for significant genes.
    logfc_col : str
        Column name with log fold change data.
    padj_col : str
        Column name with p-value data.
    gene_col : str
        Column name with gene symbol data.
    title : str | None
        PLot title, will show amount and percent of DEGs if None.
    remove_outliers : bool
        If true, remove genes with logfc value >= quartile + 15 * interquartile range.
    color_discrete_map : dict, optional
        Dict with keys "Up-regulated genes", "Down-regulated genes", "Insignificant genes"
        to assign colors.
    opacity: float
        Opacity of markers in volcano plot (0 to 1).
    return_fig : bool, optional
        Whether to return the figure object
    template : str, optional
        Plotly template name
    **kwargs : dict
        Additional arguments for px.scatter

    Returns:
    --------
    plotly.graph_objects.Figure or None
        Figure object if return_fig is True, None otherwise
    """
    template = template or pio.templates.default
    degs_df = degs_df.copy()
    degs_df = degs_df.rename(columns={logfc_col: "LFC", gene_col: "Gene", padj_col: "pval_adj"})
    if not title:
        n_degs = _calculate_n_degs(degs_df, logfc, pval)
        n_genes = len(degs_df)
        degs_percent = round(n_degs / n_genes * 100, 2)
        title = f"{n_degs} DEGs ({degs_percent}%)"
    degs_df = _process_df_for_volcano(degs_df, logfc, pval, remove_outliers=remove_outliers)
    fig = _plot_degs(
        degs_df,
        logfc,
        pval,
        title=title,
        color_discrete_map=color_discrete_map,
        opacity=opacity,
        template=template,
        **kwargs,
    )
    fig.update_xaxes(zeroline=False)
    fig.update_yaxes(zeroline=False)

    if return_fig:
        return fig
    fig.show()
    return None


def savefig(
    fig: go.Figure,
    savepath: str | Path,
    *,
    dragmode: str = "pan",
    config: dict | None = None,
    save_html: bool = True,
    save_png: bool = True,
    **kwargs,
) -> None:
    """
    Plot saving function with adjusted defaults.

    Parameters
    ----------
    fig : go.Figure
        Figure to save.
    savepath : str | Path
        Path where the figure should be saved. File extension will be added automatically
        based on the output format(s).
    dragmode : str
        Plotly dragmode setting for the figure (e.g., 'pan', 'zoom', 'select').
    config : dict | None
        Plotly config dictionary. If None, defaults to
        {"scrollZoom": True, "displaylogo": False}.
    save_html : bool
        If True, save the figure as an interactive HTML file.
    save_png : bool
        If True, save the figure as a static PNG file.
    **kwargs
        Additional keyword arguments passed to fig.write_html().

    Returns
    -------
    None
        Function saves the figure to disk and returns nothing.
    """
    savepath = Path(savepath)
    if not config:
        config = {"scrollZoom": True, "displaylogo": False}

    fig.update_layout(dragmode=dragmode)
    if save_html:
        fig.write_html(savepath.with_suffix(".html"), config=config, **kwargs)
    if save_png:
        fig.write_image(savepath.with_suffix(".png"))


import os
import platform

os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"
os.environ["MPLCONFIGDIR"] = "./tmp/matplotlib"

import scanpy as sc
import numpy as np
import pandas as pd
import yaml
import glob
import shutil
from collections import Counter

from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

# Define the path pattern for sample files
file_paths = glob.glob("tda/*.TDA_genus.txt")

# List to store individual sample DataFrames
sample_dfs = []

# Dict to store the most abundant genus of each sample
sample_genus = {}

for fp in file_paths:
    # Read the txt file
    df = pd.read_csv(fp, sep="\\t")

    # Extract sample ID from the file name
    sample_id = fp.split("/")[-1].split(".TDA_genus")[0]

    # Add a column to track the sample ID
    df["sample_id"] = sample_id

    # Append to the list
    sample_dfs.append(df)

    # Find the genus with the highest abundance
    top_genus = df.loc[df["abundance"].idxmax(), "genus"]

    # Step 4: Add to the dictionary
    sample_genus[sample_id] = top_genus

# Combine all samples into one DataFrame (long format)
combined_long = pd.concat(sample_dfs, ignore_index=True)

# top 30 genuses
genus_counts = Counter(sample_genus.values())
top30 = genus_counts.most_common(29)
top30_genus = []
for k, v in top30:
    top30_genus.append(k)

# Pivot: rows = samples, columns = genera, values = abundance (fill missing with 0)
merged_table = combined_long.pivot_table(
    index="sample_id",  # Rows: sample IDs
    columns="genus",  # Columns: genus names
    values="abundance",  # Values: abundance
    aggfunc="sum",  # Sum abundances if a genus is duplicated in the same sample
    fill_value=0,  # Fill missing genera in a sample with 0
)

adata = sc.AnnData(
    X=merged_table.values,  # Abundance matrix (samples × genera)
    obs=pd.DataFrame(index=merged_table.index),  # Sample metadata (index = sample IDs)
    var=pd.DataFrame(index=merged_table.columns),  # Genus metadata (index = genus names)
)

# Add sample metadata
df_top_genus = pd.DataFrame(  # Convert dict to DataFrame
    sample_genus.values(),  # Convert dict items to list of (key, value) tuples
    columns=["sample_genus"],  # Define column names
    index=sample_genus.keys(),
)

adata.obs = adata.obs.join(df_top_genus)  # Merge with existing sample metadata

# Perform clustering
# Build the neighborhood graph (required for both Louvain and Leiden)
n_samples = adata.n_obs
sc.tl.pca(adata, n_comps=min(n_samples // 2, 50))  # Compute PCA (n_comps ≥ n_pcs)
sc.pp.neighbors(
    adata,
    n_neighbors=int(np.sqrt(n_samples)),  # Number of nearest neighbors
    n_pcs=min(n_samples // 2, 50),  # Use PCA components
    use_rep="X_pca",  # Use PCA-reduced data (default)
)
# Run Leiden clustering
sc.tl.leiden(
    adata,
    resolution=1.0,  # Same granularity control as Louvain
    key_added="leiden",  # Store results in adata.obs["leiden"]
)

# Perform umap
sc.tl.umap(adata, random_state=0)
# Round to 10 decimal places to ensure stable hashes
adata.obsm["X_umap"] = np.round(adata.obsm["X_umap"], 10)

# Save results
adata.write_h5ad(f"umap.h5ad")
df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
df.to_pickle(f"umap.pkl")

# Plot clusters
sc.pl.umap(adata, color="sample_genus", groups=top30_genus, ncols=1, title="UMAP (Genus)", save=".pdf")
shutil.copy("figures/umap.pdf", "umap.pdf")
fig = umap(adata, color="sample_genus", template="ggplot2", ncols=1, subtitles="UMAP (Genus)", return_fig=True)
fig.write_html("figures/umap.html")
shutil.copy("figures/umap.html", "umap.html")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
