sequence = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
phe_str = rnafold(sequence);
rnaplot(phe_str, 'seq', sequence, 'format', 'tree');
title('RNA Secondary Structure Tree')