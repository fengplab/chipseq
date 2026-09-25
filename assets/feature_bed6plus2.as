table featureBed6Plus2
"RepeatMasker / CenSat features prepared by nf-core/chipseq (features_to_bed.py)"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Feature name (e.g. repeat name or CenSat annotation)"
    uint   score;      "Score (0-1000)"
    char[1] strand;    "+, - or ."
    string class;      "Feature class (e.g. RepeatMasker repClass or CenSat satellite type)"
    string family;     "Feature family (e.g. RepeatMasker repFamily)"
    )
