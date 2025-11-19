#include "x3dna.h"
typedef struct {
    char inpfile[BUF512];
    char outfile[BUF512];
    double scale_factor;
} struct_args;
static void pdb2img_usage(void)
{
    help3dna_usage("pdb2img");
}

static void set_defaults(struct_args *args)
{
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "stdout");
    args->scale_factor = 0.0;
}

static void pdb2img_cmdline(int argc, char *argv[], struct_args *args, long *opts)
{
    long i, j;
    if (argc < 3)
        pdb2img_usage();
    set_defaults(args);
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        if (check_global_options(argv[i]))
            continue;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)
            if (argv[i][j] == 'F')
                opts[0] = 1;
            else if (argv[i][j] == 'R')
                opts[0] = 2;
            else if (argv[i][j] == 'B')
                opts[1] = 1;
            else if (argv[i][j] == 'C')
                opts[2] = 1;
            else if (argv[i][j] == 'P')
                opts[4] = 1;
            else if (argv[i][j] == 'I')
                opts[7] = 1;
            else if (argv[i][j] == 'N')
                opts[8] = 1;
            else if (argv[i][j] == 'U')
                opts[9] = 1;
            else if (argv[i][j] == 'M')
                opts[10] = 1;
            else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", &args->scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    pdb2img_usage();
                } else
                    break;
            } else
                pdb2img_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else if (argc == i + 1)
        strcpy(args->inpfile, argv[i]);
    else
        pdb2img_usage();
    if (opts[10])
        opts[7] = 1;
}

static void pdb2alc(char *inpfile, char *alcfile, long peptide)
{
    char BDIR[BUF512];
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **orien, **org, **xyz;
    long np, num, num_residue;
    long *cidx, *mchain, *resid, *ResSeq, *RY, **seidx;
    num = number_of_atoms(inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 4);
    if (peptide) {
        get_BDIR(BDIR, SNAP_PEP_PDB);
        cidx = lvector(1, num_residue);
        resid = lvector(1, num_residue);
        np = 4 * num_residue;
        mchain = lvector(1, np);
        peptide_info(num_residue, seidx, AtomName, ResName, ChainID,
                     ResSeq, Miscs, xyz, resid, cidx, mchain);
        peptide_frame(num_residue, BDIR, resid, mchain, xyz, orien, org);
        peptide_blks(num_residue, BDIR, cidx, orien, org, alcfile);
        free_lvector(cidx, 1, num_residue);
        free_lvector(resid, 1, num_residue);
        free_lvector(mchain, 1, np);
    } else {
        bseq = cvector(1, num_residue);
        RY = lvector(1, num_residue);
        get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);
        get_BDIR(BDIR, "Atomic_A.pdb");
        base_frame(num_residue, bseq, seidx, RY, AtomName, ResName,
                   ChainID, ResSeq, Miscs, xyz, BDIR, orien, org);
        get_BDIR(BDIR, "Block_R.alc");
        base_blks(num_residue, RY, orien, org, bseq, BDIR, alcfile);
        free_cvector(bseq, 1, num_residue);
        free_lvector(RY, 1, num_residue);
    }
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 4);
}

int main(int argc, char *argv[])
{
    char alcfile[BUF512];
    struct_args args;
    long opts[BUF512] = { 0 };
    remove_file(BBLKALC_FILE);
    remove_file(PBLKALC_FILE);
    set_my_globals(argv[0]);
    pdb2img_cmdline(argc, argv, &args, opts);
    (opts[4]) ? strcpy(alcfile, PBLKALC_FILE) : strcpy(alcfile, BBLKALC_FILE);
    pdb2alc(args.inpfile, alcfile, opts[4]);
    process_alc(alcfile, args.outfile, args.scale_factor, opts);
    clear_my_globals();
    return 0;
}
