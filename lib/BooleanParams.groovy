class BooleanParams {
    static final List<String> NAMES = [
        'help', 'minimeta', 'prepare_databases', 'monochrome_logs',
        'igenomes_ignore', 'show_hidden_params', 'validationShowHiddenParams',
        'validationSchemaIgnoreParams', 'plaintext_email',
        'single_end', 'bulk', 'mg', 'notrim', 'saveTrimmed',
        'saveAlignedIntermediates', 'allow_multi_align', 'ass',
        'no_normalize', 'euk', 'fungus', 'saturation', 'snv', 'cnv',
        'doubletd', 'bbmap', 'kraken', 'genomad', 'checkm2', 'blastn',
        'blob', 'acdc', 'pangenome', 'completeness', 'tree', 'eggnog',
        'kofam', 'acquired', 'point', 'split', 'split_euk', 'graphbin',
        'gtdbtk', 'pasa', 'run_cooccurrence_checkm', 'nanopore', 'remap'
    ].asImmutable()

    static boolean value(Map params, String name) {
        def raw = params[name]
        if (raw instanceof Boolean) return raw
        if (raw instanceof CharSequence && ['true', 'false'].contains(raw.toString().toLowerCase())) {
            return raw.toString().toBoolean()
        }
        throw new IllegalArgumentException("--${name} must be true or false (got: ${raw})")
    }

    static void validate(Map params) {
        NAMES.each { name ->
            value(params, name)
        }
    }
}
