## Alignment reorder

You may wish to reorder the samples in the alignment, such as to group similar samples together to view haplotypes in a program such as MEGA or aliview.

First, create a list of the order you wish for the samples to be in. If you have created a tree, run this command:

```bash
cat $output_tree.treefile | grep -oP '[^():,;]+' | grep -v '^[0-9.]\+$' > sample_order.txt
```

Afterwards, use this code to reorder the fasta sequences:

```bash
awk -v order="sample_order.txt" '
BEGIN {
    while ((getline line < order) > 0) {
        seq_order[line] = ""
        seq_names[++n] = line
    }
    close(order)
}
/^>/ {
    if(seq_name) fasta[seq_name]=seq_seq
    seq_name=substr($0,2)
    seq_seq=""
    next
}
{
    seq_seq=seq_seq $0
}
END {
    if(seq_name) fasta[seq_name]=seq_seq
    for(i=1;i<=n;i++) {
        name=seq_names[i]
        if(name in fasta) {
            print ">" name
            print fasta[name]
        }
    }
}
' alignment.fa > reordered.fa
```
