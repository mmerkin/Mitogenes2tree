#!/usr/bin/awk -f
# Usage: ./update_genbank.awk GENUS SPECIES input.gb > updated.gb

BEGIN {
    if (ARGC < 4) {
        print "usage: ./update_genbank.awk GENUS SPECIES input.gb > updated.gb" > "/dev/stderr"
        exit 1
	}
    GENUS = ARGV[1]
    SPECIES = ARGV[2]
    ARGV[1] = ""; ARGV[2] = ""

    genome_length = 0
    last_feature_end = 0

    feature_indent = "     "
    qualifier_indent = "                     "

    current_feature = ""
    last_gene_name = ""

    in_tRNA = 0
    tRNA_has_gene = 0

    pending_rrna = 0
    rrna_line = ""
    rrna_location = ""
}

# Replace Locus

/^LOCUS/ {
    match($2, /[0-9]+/, arr)
    genome_length = arr[0]

    printf "%s       %s       %s bp    DNA     circular     %s\n",
        $1, GENUS"_"SPECIES"_mtDNA", genome_length, $6
    next
}


# Replace DEFINITION line
/^DEFINITION/ {
    print "DEFINITION  " GENUS " " SPECIES " mitochondrion, complete genome."
    next
}

# Replace Unknown. in SOURCE
/^[[:space:]]*SOURCE[[:space:]]+Unknown\./ { sub(/Unknown\./, "mitochondrion " GENUS " " SPECIES) }

# Replace Unknown. in ORGANISM
/^  ORGANISM/ { sub(/Unknown\./, GENUS " " SPECIES) }

# Update source feature

/^\s*\/mol_type="genomic DNA"/ {
    print "                     /organism=\"" GENUS " " SPECIES "\""
    print "                     /mol_type=\"mitochondrial DNA\""
    next
}

# Detect feature lines
/^     [A-Za-z_-]+[[:space:]]+/ {

    # Close tRNA block if needed
    if (in_tRNA && !tRNA_has_gene && last_gene_name ~ /^trn/) {
        print qualifier_indent "/gene=\"" last_gene_name "\""
    }

    in_tRNA = 0
    tRNA_has_gene = 0

    current_feature = $1

    if (match($0, /([0-9]+)\.\.([0-9]+)/, arr))
        last_feature_end = arr[2]

    # Buffer rRNA
    if (current_feature == "rRNA") {
        pending_rrna = 1
        rrna_line = $0
        rrna_location = $0
        sub(/^[[:space:]]*rRNA[[:space:]]+/, "", rrna_location)
        next
    }

    print
    if (current_feature == "tRNA")
        in_tRNA = 1

    next
}

# Capture gene name
current_feature == "gene" && /\/gene="/ {
    match($0, /"([^"]+)"/, arr)
    last_gene_name = arr[1]
    print
    next
}

# FIX: Correct wrong rrn gene inside tRNA
in_tRNA && /\/gene="rrn[LS]"/ {
    print qualifier_indent "/gene=\"" last_gene_name "\""
    tRNA_has_gene = 1
    next
}

# Detect correct gene inside tRNA
in_tRNA && /\/gene="/ {
    tRNA_has_gene = 1
    print
    next
}

# Repair rRNA block
pending_rrna && /\/product="/ {

    rrna_gene = ""
    if ($0 ~ /l-rRNA/) rrna_gene = "rrnL"
    else if ($0 ~ /s-rRNA/) rrna_gene = "rrnS"

    print feature_indent "gene            " rrna_location
    print qualifier_indent "/gene=\"" rrna_gene "\""

    print rrna_line
    print qualifier_indent "/gene=\"" rrna_gene "\""
    print

    pending_rrna = 0
    next
}

# Insert D-loop
/^BASE COUNT/ {

    if (in_tRNA && !tRNA_has_gene && last_gene_name ~ /^trn/) {
        print qualifier_indent "/gene=\"" last_gene_name "\""
    }

    if (genome_length > 0 && last_feature_end > 0) {
        print feature_indent "D-Loop          " (last_feature_end+1) ".." genome_length
    }

    print
    next
}

# Default
{ print }
