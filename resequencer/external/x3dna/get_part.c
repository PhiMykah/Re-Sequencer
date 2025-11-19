#include "x3dna.h"
typedef struct {
    char inpfile[BUF512];
    char outfile[BUF512];
    long na;
    long protein;
    long ligand;
    long dwater;
    long dhatom;
    long bb;
    long msc;
    long all;
} struct_args;
static void getpart_usage(void)
{
    help3dna_usage("get_part");
}

static void set_defaults(struct_args *args)
{
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->na = FALSE;
    args->protein = FALSE;
    args->ligand = FALSE;
    args->dwater = FALSE;
    args->dhatom = FALSE;
    args->bb = 0;
    args->msc = FALSE;
    args->all = FALSE;
}

static void getpart_cmdline(int argc, char *argv[], struct_args *args)
{
    long i, j;
    if (argc < 3)
        getpart_usage();
    set_defaults(args);
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        if (check_global_options(argv[i]))
            continue;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)
            if (argv[i][j] == 'N')
                args->na = TRUE;
            else if (argv[i][j] == 'P')
                args->protein = TRUE;
            else if (argv[i][j] == 'C')
                args->msc = TRUE;
            else if (argv[i][j] == 'T')
                args->ligand = TRUE;
            else if (argv[i][j] == 'W')
                args->dwater = TRUE;
            else if (argv[i][j] == 'D')
                args->dhatom = TRUE;
            else if (argv[i][j] == 'L')
                args->bb = 1;
            else if (argv[i][j] == 'B')
                args->bb = 2;
            else if (argv[i][j] == 'X')
                args->bb = 3;
            else if (argv[i][j] == 'Z')
                args->bb = 4;
            else if (argv[i][j] == 'A')
                args->all = TRUE;
            else
                getpart_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        getpart_usage();
    if (!args->na && !args->protein && !args->ligand && !args->bb && !args->dwater
        && !args->dhatom)
        args->na = TRUE;
    if (args->bb)
        args->na = FALSE;
}

static long is_ringatom(char *aname, long resid)
{
    static char *RingAtom[] = { RA_LIST };
    long i, rnum;
    rnum = (resid) ? 9 : 6;
    for (i = 0; i < rnum; i++)
        if (!strcmp(aname, RingAtom[i]))
            return 1;
    return 0;
}

static void check_color_name(char *chnname, char *colname, long mscol_num, char **mscol_name)
{
    long k;
    if ((!strcmp(chnname, "protein") || !strcmp(chnname, "nucleic_acid")) &&
        !strcmp(colname, "by_chain"))
        return;
    if (strcmp(chnname, "ligand") && strcmp(chnname, "protein") &&
        strcmp(chnname, "nucleic_acid") && strlen(chnname) != 1)
        fatal("Wrong color format: %s - %s\n", chnname, colname);
    k = num_strmatch(colname, mscol_name, 1, mscol_num);
    if (!k)
        fatal("Color name <%s> is NOT legal in MolScript\n", colname);
    if (k > 1)
        fprintf(stderr, "Color name <%s> is repeated in MolScript\n", colname);
}

static long molscript_colors(char **mscol)
{
    char *COL_MolScript = "col_mname.dat";
    char BDIR[BUF512], str[BUF512], tmp[BUF512];
    long mc_num = 0;
    FILE *fp;
    get_BDIR(BDIR, COL_MolScript);
    strcat(BDIR, COL_MolScript);
    fprintf(stderr, " ...... reading file: %s ...... \n", COL_MolScript);
    fp = open_file(BDIR, "r");
    while ((fgets(str, sizeof str, fp) != NULL)) {
        if (str[0] == '#')
            continue;
        lowerstr(str);
        if (sscanf(str, "%s", tmp) != 1)
            fatal("can not read color name\n");
        if (strlen(tmp) > MAXCLEN)
            fatal("color name too long\n");
        if (++mc_num > BUF512)
            fatal("too many colors\n");
        strcpy(mscol[mc_num], tmp);
    }
    close_file(fp);
    return mc_num;
}

static long chain_colors(char *chains, char **ccolors, char **pnl, double *NA_coilradius)
{
    char *COL_CHAIN_FILE = "col_chain.dat";
    char BDIR[BUF512], str[BUF512], chnname[BUF512], colname[BUF512];
    char **mscol_name;
    long mscol_num, chain_num = 0;
    FILE *fp;
    mscol_name = cmatrix(1, BUF512, 0, MAXCLEN);
    mscol_num = molscript_colors(mscol_name);
    get_BDIR(BDIR, COL_CHAIN_FILE);
    strcat(BDIR, COL_CHAIN_FILE);
    fprintf(stderr, " ...... reading file: %s ...... \n", COL_CHAIN_FILE);
    chains[0] = '*';
    strcpy(ccolors[0], "goldenrod");
    strcpy(pnl[0], "by_chain");
    strcpy(pnl[1], "by_chain");
    fp = open_file(BDIR, "r");
    while ((fgets(str, sizeof str, fp) != NULL)) {
        if (str[0] == '#')
            continue;
        lowerstr(str);
        if (sscanf(str, "%s %s", chnname, colname) != 2)
            continue;
        if (chnname[0] == '#' || colname[0] == '#')
            continue;
        check_color_name(chnname, colname, mscol_num, mscol_name);
        if (!strcmp(chnname, "protein"))
            strcpy(pnl[0], colname);
        else if (!strcmp(chnname, "nucleic_acid")) {
            strcpy(pnl[1], colname);
            if (sscanf(str, "%*s %*s %lf", NA_coilradius) == 1) {
                fprintf(stderr, "read in coil radius: %.2f\n", *NA_coilradius);
                if (*NA_coilradius < 0) {
                    *NA_coilradius = -*NA_coilradius;
                    fprintf(stderr, "          changed to: %.2f\n", *NA_coilradius);
                }
            } else
                fprintf(stderr, "use default coil radius: %.2f\n", *NA_coilradius);
        }
        if (strlen(chnname) > 1)
            continue;
        if (chnname[0] == '*') {
            chains[0] = chnname[0];
            strcpy(ccolors[0], colname);
            continue;
        }
        if (++chain_num >= BUF512)
            fatal("too many chain ID - color pairs\n");
        chains[chain_num] = toupper((int) chnname[0]);
        strcpy(ccolors[chain_num], colname);
    }
    close_file(fp);
    chains[chain_num + 1] = '\0';
    free_cmatrix(mscol_name, 1, BUF512, 0, MAXCLEN);
    return chain_num;
}

static long expt_chainID(long num_residue, long **seidx, long *res_type, char *ChainID,
                         char *Echains, char *chain_type)
{
    char c1;
    long i, numc = 0;
    Echains[0] = '\0';
    for (i = 1; i <= num_residue; i++) {
        c1 = ChainID[seidx[i][1]];
        if (c1 == ' ')
            continue;
        if (strchr(Echains, c1) == NULL) {
            Echains[numc] = c1;
            if (res_type[i] >= 0)
                chain_type[numc] = 'N';
            else if (res_type[i] == -1)
                chain_type[numc] = 'P';
            else if (res_type[i] == -3)
                chain_type[numc] = 'L';
            else
                chain_type[numc] = 'O';
            if (++numc > BUF512)
                fatal("too many chain IDs in this PDB file\n");
            Echains[numc] = '\0';
        }
    }
    chain_type[numc] = '\0';
    fprintf(stderr, "Expt. chains (%ld)\n", numc);
    fprintf(stderr, "             %s\n", Echains);
    fprintf(stderr, "             %s\n", chain_type);
    return numc;
}

static void molauto_protein(char *chains, char **ccolors, char *Echains, char *chain_type,
                            char **pnl, FILE *fp)
{
    char *molauto_file = "temp";
    char c1, c0 = '\0', *pchar, str[BUF512], ccol[BUF512];
    long j, k;
    FILE *fpmol;
    fprintf(fp, "! protein part\n");
    fpmol = open_file(molauto_file, "r");
    while ((fgets(str, sizeof str, fpmol) != NULL)) {
        pchar = strstr(str, " from ");
        if (pchar == NULL || strstr(str, " to ") == NULL)
            continue;
        k = pchar - str + 6;
        c1 = str[k];
        pchar = strchr(Echains, c1);
        if (pchar == NULL) {
            fprintf(stderr, "Chain ID %c <%s> does NOT exists in PDB file\n", c1, molauto_file);
            j = 0;
        } else {
            k = pchar - Echains;
            if (chain_type[k] != 'P')
                fprintf(stderr, "Chain ID %c <%s> is NOT a protein chain\n", c1, molauto_file);
            pchar = strchr(chains, c1);
            j = (pchar == NULL) ? 0 : pchar - chains;
        }
        if (c1 != c0) {
            strcpy(ccol, (!strcmp(pnl[0], "by_chain")) ? ccolors[j] : pnl[0]);
            fprintf(fp, "\n    set planecolour %s;\n", ccol);
            c0 = c1;
        }
        fprintf(fp, "  %s", str);
    }
    close_file(fpmol);
}

static void molscript_na(char *chains, char **ccolors, long num_Echains, char *Echains,
                         char *chain_type, char **pnl, double NA_coilradius, FILE *fp)
{
    char *pchar, ccol[BUF512];
    long i, j;
    fprintf(fp, "\n! Nucleic acid colored by chain IDs\n\n");
    fprintf(fp, "    set coilradius %.2f;\n", NA_coilradius);
    for (i = 0; i < num_Echains; i++) {
        if (chain_type[i] != 'N')
            continue;
        pchar = strchr(chains, Echains[i]);
        j = (pchar == NULL) ? 0 : pchar - chains;
        strcpy(ccol, (!strcmp(pnl[1], "by_chain")) ? ccolors[j] : pnl[1]);
        fprintf(fp, "    set planecolour %s;\n", ccol);
        fprintf(fp, "    double-helix chain %c;\n", Echains[i]);
    }
    if (!num_Echains) {
        strcpy(ccol, (!strcmp(pnl[1], "by_chain")) ? ccolors[0] : pnl[1]);
        fprintf(fp, "    set planecolour %s;\n", ccol);
        fprintf(fp, "    double-helix nucleotides;\n");
    }
}

static void msc_setting(char *inpfile, char *outfile, long *res_type, char *ChainID,
                        long num_residue, long **seidx)
{
    char chains[BUF512 + 1], Echains[BUF512 + 1], **ccolors, **pnl;
    char chain_type[BUF512 + 1];
    double NA_coilradius = 0.5;
    long i, num_Echains;
    FILE *fp;
    ccolors = cmatrix(0, BUF512, 0, MAXCLEN);
    pnl = cmatrix(0, 1, 0, MAXCLEN);
    (void) chain_colors(chains, ccolors, pnl, &NA_coilradius);
    num_Echains = expt_chainID(num_residue, seidx, res_type, ChainID, Echains, chain_type);
    fp = open_file(outfile, "w");
    fprintf(fp, "! MolScript input file generated by 3DNA\n\n");
    fprintf(fp, "plot\n\n");
    fprintf(fp, "    read mol \"%s\";\n", inpfile);
    fprintf(fp, "    set colourparts off;\n");
    fprintf(fp, "    set segments 6, plane2colour grey 0.8;\n\n");
    molauto_protein(chains, ccolors, Echains, chain_type, pnl, fp);
    for (i = 1; i <= num_residue; i++)
        if (res_type[i] >= 0) {
            molscript_na(chains, ccolors, num_Echains, Echains, chain_type, pnl,
                         NA_coilradius, fp);
            break;
        }
    fprintf(fp, "\nend_plot\n");
    close_file(fp);
    free_cmatrix(ccolors, 0, BUF512, 0, MAXCLEN);
    free_cmatrix(pnl, 0, 1, 0, MAXCLEN);
}

static void write_selected(char *inpfile, char *outfile, long num_residue, long **seidx,
                           long *res_type, char **AtomName, char **ResName, char *ChainID,
                           long *ResSeq, double **xyz, char **Miscs, struct_args *args)
{
    char *str;
    long i, iok, j, k, ib, ie, nidx = 0, resid;
    long outnum = 0;
    FILE *fp;
    fp = open_file(outfile, "w");
    print_pdb_title(inpfile, "*", fp);
    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        resid = res_type[i];
        if (args->bb && resid >= 0)
            for (j = ib; j <= ie; j++) {
                iok = 0;
                k = is_baseatom(AtomName[j]);
                if (args->bb == 1) {
                    if (!k)
                        iok = 1;
                } else if (args->bb == 2) {
                    if (k)
                        iok = 1;
                } else {
                    if (k) {
                        if (is_ringatom(AtomName[j], resid))
                            iok = 1;
                    } else {
                        if (args->bb == 3 && (strcmp(AtomName[j], " O1P") &&
                                              strcmp(AtomName[j], " O2P")))
                            iok = 1;
                        if (args->bb == 4 && !strcmp(AtomName[j], " P  "))
                            iok = 1;
                    }
                }
                if (!iok)
                    continue;
                str = Miscs[j];
                fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                        (str[0] == 'A') ? "ATOM  " : "HETATM", ++nidx, AtomName[j],
                        str[1], ResName[j], ChainID[j], ResSeq[j], str[2], xyz[j][1],
                        xyz[j][2], xyz[j][3], str + 3);
                outnum++;
            }
        k = 0;
        if (args->dwater && resid != -6)
            k = 1;
        if (args->na && resid >= 0)
            k = 1;
        if (args->protein && resid == -1)
            k = 1;
        if (args->ligand && resid == -3)
            k = 1;
        if (k) {
            pdb_record(ib, ie, &nidx, 0, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, fp);
            outnum++;
        }
    }
    fprintf(fp, "END\n");
    close_file(fp);
    if (args->dhatom) {
        cpcat_file(outnum ? outfile : inpfile, TMP_FILE, "copy");
        delH_pdbfile(TMP_FILE, outfile);
    }
}

int main(int argc, char *argv[])
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    double **xyz;
    long num, num_residue;
    long *ResSeq, **seidx, *res_type;
    struct_args args;
    set_my_globals(argv[0]);
    getpart_cmdline(argc, argv, &args);
    num = number_of_atoms(args.inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(args.inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    if (args.all) {
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, args.outfile);
        free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
        return 0;
    }
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    res_type = lvector(1, num_residue);
    residue_wtype(num_residue, seidx, ResName, AtomName, xyz, Miscs, res_type, FALSE);
    if (args.msc)
        msc_setting(args.inpfile, args.outfile, res_type, ChainID, num_residue, seidx);
    else
        write_selected(args.inpfile, args.outfile, num_residue, seidx, res_type, AtomName,
                       ResName, ChainID, ResSeq, xyz, Miscs, &args);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lvector(res_type, 1, num_residue);
    clear_my_globals();
    return 0;
}
