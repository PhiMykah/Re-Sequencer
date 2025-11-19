#include "x3dna.h"
typedef struct {
    char hbfile[BUF512];
    char lkgfile[BUF512];
    char inpfile[BUF512];
    char outfile[BUF512];
    double scale_factor;
    double ballrad;
} struct_args;
static void stack2img_usage(void)
{
    help3dna_usage("stack2img");
}

static void set_defaults(struct_args *args)
{
    strcpy(args->hbfile, "");
    strcpy(args->lkgfile, "");
    strcpy(args->inpfile, "");
    strcpy(args->outfile, "");
    args->scale_factor = 0.0;
    args->ballrad = 0.0;
}

static void stack2img_cmdline(int argc, char *argv[], struct_args *args, long *opts)
{
    char str[BUF512];
    long i, j;
    if (argc < 3)
        stack2img_usage();
    set_defaults(args);
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        if (check_global_options(argv[i]))
            continue;
        strcpy(str, argv[i]);
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)
            if (argv[i][j] == 'F')
                opts[0] = 1;
            else if (argv[i][j] == 'R')
                opts[0] = 2;
            else if (argv[i][j] == 'D')
                opts[1] = 1;
            else if (argv[i][j] == 'O')
                opts[2] = 1;
            else if (argv[i][j] == 'L')
                opts[3] = 1;
            else if (argv[i][j] == 'C')
                opts[4] = 1;
            else if (argv[i][j] == 'B')
                opts[5] = 1;
            else if (argv[i][j] == 'T')
                opts[6] = 1;
            else if (argv[i][j] == 'M')
                opts[7] = 1;
            else if (argv[i][j] == 'P')
                opts[8] = 1;
            else if (argv[i][j] == 'X') {
                if (sscanf(str, "-X=%s", args->hbfile) != 1 &&
                    sscanf(str, "-x=%s", args->hbfile) != 1) {
                    fprintf(stderr, "wrong format for reading H-bond list\n");
                    stack2img_usage();
                } else {
                    opts[1] = 1;
                    break;
                }
            } else if (argv[i][j] == 'Y') {
                if (sscanf(str, "-Y=%s", args->lkgfile) != 1 &&
                    sscanf(str, "-y=%s", args->lkgfile) != 1) {
                    fprintf(stderr, "wrong format for reading bond list\n");
                    stack2img_usage();
                } else
                    break;
            } else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", &args->scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    stack2img_usage();
                } else
                    break;
            } else if (argv[i][j] == 'V') {
                if (sscanf(argv[i], "-V=%lf", &args->ballrad) != 1) {
                    fprintf(stderr, "wrong format for atom radius setting\n");
                    stack2img_usage();
                } else
                    break;
            } else
                stack2img_usage();
    }
    if (argc == i + 2) {
        strcpy(args->inpfile, argv[i]);
        strcpy(args->outfile, argv[i + 1]);
    } else
        stack2img_usage();
    if (opts[6] && opts[7]) {
        fprintf(stderr, "Both -t and -m supplied: ignoring -m\n");
        opts[7] = 0;
    }
}

static void get_pair_num(char *inpfile, long **pair_num)
{
    char str[BUF512];
    char *item[MBASES];
    FILE *fp;
    long i, j, k, nitem;
    fp = open_file(inpfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        nitem = upperstr(str);
        if (strstr(str, "REMARK") && (strstr(str, "LOWER:") || strstr(str, "UPPER:"))) {
            j = strstr(str, "LOWER:") ? 1 : 2;
            nitem = itemize(str, item, MBASES);
            for (i = 2; i <= nitem; i++) {
                if (sscanf(item[i], "%ld", &k) != 1)
                    break;
                if (k <= 0)
                    fatal("residue number out of range (%ld <= 0)!\n", k);
                pair_num[j][i - 2] = k;
            }
            if (i <= nitem)
                fatal("containing too few base residue numbers\n");
            if (j == 2)
                break;
        }
    }
    close_file(fp);
}

static void stack2fig(long nobj, long *idx, long *aidx, long *depth, long **allobj,
                      double **xyz, long **linkage, long **hb_linkage, long **ring_atom,
                      long *ibase, long is_color, long **pair_num, long label_ring,
                      long *ResSeq, double ballrad, FILE *fp)
{
    static char *cmn_base = CX_LIST;
    static char *label_font = "0 18";
    static char *format = "%6.0f%6.0f";
    static long FTSIZE = 18;
    static double FTFACTOR = 0.4;
    double p2f = 1200.0 / 72.0;
    double dot_sep, msat, Msat, xc, yc, yoffset, vdw_radii[NELE];
    long dlcol, dwidth, line_width, join_style, cap_style, o_sides, mfcol;
    long bcol_code = 0, i, i1, i2, ib = 0, j, k, Mfill, mfill, rnum, bpidx;
    long bpwidth[3], **bc_idx;
    atom_info(3, NULL, NULL, vdw_radii);
    if (ballrad > 0.0) {
        fprintf(fp, "0 100 #FF1493\n");
        fprintf(fp, "0 101 #C8C8C8\n");
        fprintf(fp, "0 102 #F00000\n");
        fprintf(fp, "0 103 #404040\n");
        fprintf(fp, "0 104 #8F8FFF\n");
        fprintf(fp, "0 105 #FFC832\n");
        fprintf(fp, "0 106 #FFA500\n");
        fprintf(fp, "0 107 #DAA520\n");
        fprintf(fp, "0 108 #00FF00\n");
        fprintf(fp, "0 109 #A52A2A\n");
        fprintf(fp, "0 110 #A020F0\n");
        fprintf(fp, "0 111 #DAA520\n\n");
    }
    bc_idx = lmatrix(0, 1, 0, 6);
    get_fig_pars(&dot_sep, &dlcol, &dwidth, &bpwidth[1], &bpwidth[2], bc_idx,
                 &msat, &Msat, &o_sides, &line_width, &join_style, &cap_style, &mfcol);
    yoffset = FTFACTOR * FTSIZE * p2f;
    for (i = 1; i <= nobj; i++) {
        j = idx[i];
        if (allobj[1][j] == 0) {
            ib = allobj[2][j];
            k = (is_color) ? 100 + aidx[ib] : 0;
            mfill = (is_color) ? 20 : 15;
            fprintf(fp, "1 3 0 %ld %ld %ld %ld 0 %ld 0.0 1 0.0 ", bpwidth[1], k, k,
                    depth[i], mfill);
            fprintf(fp, format, xyz[ib][1], xyz[ib][2]);
            xc = ballrad * vdw_radii[aidx[ib]] * p2f;
            fprintf(fp, " %.0f %.0f 0 0 0 0\n", xc, xc);
            continue;
        }
        if (allobj[1][j] == -2)
            bpidx = (allobj[2][j] < 0) ? 2 : 1;
        else {
            ib = allobj[2][j];
            bcol_code = bc_idx[is_color][ibase[ib]];
            bpidx = (lval_in_set(ib, 1, pair_num[2][0], pair_num[2])) ? 2 : 1;
        }
        if (allobj[1][j] > 0) {
            fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0 %2ld %2ld 0 0 0 2\n",
                    bpwidth[bpidx], bcol_code, depth[i], join_style, cap_style);
            for (k = 1; k <= 2; k++) {
                i1 = linkage[allobj[1][j]][k];
                fprintf(fp, format, xyz[i1][1], xyz[i1][2]);
            }
            fprintf(fp, "\n");
        } else if (allobj[1][j] == -1) {
            k = lround(msat * 20);
            mfill = (bcol_code) ? 20 + k : 20 - k;
            k = lround(Msat * 20);
            Mfill = (bcol_code) ? 40 - k : k;
            i1 = allobj[2][j];
            rnum = ring_atom[i1][10];
            fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 ",
                    bpwidth[bpidx], bcol_code, bcol_code, depth[i]);
            if (bpidx == 2)
                fprintf(fp, "%2ld", mfill);
            else
                fprintf(fp, "%2ld", Mfill);
            fprintf(fp, " 0.0 %2ld %2ld 0 0 0 %2ld\n", join_style, cap_style, rnum + 1);
            for (k = 1; k <= rnum; k++) {
                i2 = ring_atom[i1][k];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                if (!(k % 4))
                    fprintf(fp, "\n");
            }
            i2 = ring_atom[i1][1];
            fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
            fprintf(fp, "\n");
            if (rnum == 9) {
                fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld", bpwidth[bpidx], bcol_code, depth[i]);
                fprintf(fp, " 0 -1 0.0 %2ld %2ld 0 0 0 2\n", join_style, cap_style);
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                i2 = ring_atom[i1][6];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                fprintf(fp, "\n");
            }
            if (label_ring) {
                fprintf(fp, "4 1 0 %4ld %s %ld 0.0 4 165 165", depth[i] - 1, label_font,
                        FTSIZE);
                if (Gvars.label_RC8_YC6) {
                    k = (rnum == 9) ? ring_atom[i1][8] : ring_atom[i1][5];
                    xc = xyz[k][1];
                    yc = xyz[k][2];
                } else {
                    xc = 0.5 * (xyz[ring_atom[i1][1]][1] + xyz[ring_atom[i1][4]][1]);
                    yc = 0.5 * (xyz[ring_atom[i1][1]][2] + xyz[ring_atom[i1][4]][2]);
                }
                fprintf(fp, format, xc, yc + yoffset);
                fprintf(fp, " %c%ld\\001\n", cmn_base[ibase[ib]], ResSeq[i2]);
            }
        } else {
            if (!is_color)
                dlcol = 0;
            fprintf(fp, "2 1 2 %2ld %2ld 0", bpwidth[bpidx], dlcol);
            fprintf(fp, " %4ld 0 -1 %5.1f %2ld %2ld 0 0 0 2\n",
                    depth[i], dot_sep, join_style, cap_style);
            for (k = 2; k <= 3; k++) {
                i2 = hb_linkage[labs(allobj[2][j])][k];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
            }
            fprintf(fp, "\n");
        }
    }
    free_lmatrix(bc_idx, 0, 1, 0, 6);
}

static void stack2ps(long nobj, long *idx, long *aidx, long **allobj, double **xyz,
                     long **linkage, long **hb_linkage, long **ring_atom, long *ibase,
                     long is_color, long **pair_num, long label_ring, long *ResSeq,
                     double ballrad, FILE *fp_ok)
{
    char str[BUF512], bname = ' ', *format = "%7.1f%7.1f";
    char *temp_eps = "_eps_temp.txt", *temp_labels = "_labels_temp.txt";
    static char bc_idx[2][8] = { "XXXXXXX",
        CX_LIST
    };
    double vdw_radii[NELE], xc, yc;
    long bpidx, i, i1, i2 = 0, ib = 0, j, k, rnum;
    FILE *fp;
    atom_info(3, NULL, NULL, vdw_radii);
    fp = open_file(temp_eps, "w");
    for (i = 1; i <= nobj; i++) {
        j = idx[i];
        if (allobj[1][j] == 0) {
            ib = allobj[2][j];
            fprintf(fp, "NP ");
            fprintf(fp, format, xyz[ib][1], xyz[ib][2]);
            fprintf(fp, " %.2f CIRCLE ", ballrad * vdw_radii[aidx[ib]]);
            if (is_color)
                fprintf(fp, "ATOM%2.2ld setrgbcolor fill\n", aidx[ib]);
            else
                fprintf(fp, "GREYLEVEL setgray fill\n");
            continue;
        }
        if (allobj[1][j] == -2)
            bpidx = (allobj[2][j] < 0) ? 2 : 1;
        else {
            ib = allobj[2][j];
            bname = bc_idx[is_color][ibase[ib]];
            bpidx = (lval_in_set(ib, 1, pair_num[2][0], pair_num[2])) ? 2 : 1;
        }
        (bpidx == 2) ? fprintf(fp, "NP W2 ") : fprintf(fp, "NP W1 ");
        if (allobj[1][j] > 0) {
            fprintf(fp, "%cl ", bname);
            for (k = 1; k <= 2; k++) {
                i1 = linkage[allobj[1][j]][k];
                fprintf(fp, format, xyz[i1][1], xyz[i1][2]);
            }
            fprintf(fp, " LN\n");
        } else if (allobj[1][j] == -1) {
            fprintf(fp, "%cl\n", bname);
            i1 = allobj[2][j];
            rnum = ring_atom[i1][10];
            for (k = 1; k <= rnum; k++) {
                i2 = ring_atom[i1][k];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                if (!(k % 4))
                    fprintf(fp, "\n");
            }
            fprintf(fp, " R%1ld\n", rnum);
            if (bpidx == 2)
                fprintf(fp, "  gsave %cm grestore stroke\n", bname);
            else
                fprintf(fp, "  gsave %cM grestore stroke\n", bname);
            if (rnum == 9) {
                fprintf(fp, "NP ");
                i2 = ring_atom[i1][1];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                i2 = ring_atom[i1][6];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
                fprintf(fp, " LN\n");
            }
            if (label_ring) {
                if (is_color)
                    fprintf(fp, "0.618 setgray ");
                else {
                    fprintf(fp, "0.918 setgray ");
                    bname = bc_idx[1][ibase[ib]];
                }
                if (Gvars.label_RC8_YC6) {
                    k = (rnum == 9) ? ring_atom[i1][8] : ring_atom[i1][5];
                    xc = xyz[k][1];
                    yc = xyz[k][2];
                } else {
                    xc = 0.5 * (xyz[ring_atom[i1][1]][1] + xyz[ring_atom[i1][4]][1]);
                    yc = 0.5 * (xyz[ring_atom[i1][1]][2] + xyz[ring_atom[i1][4]][2]);
                }
                fprintf(fp, format, xc, yc);
                fprintf(fp, " YOFFSET sub moveto (%c%ld) SCENTER\n", bname, ResSeq[i2]);
            }
        } else {
            (is_color) ? fprintf(fp, "Dl ") : fprintf(fp, "Xl ");
            for (k = 2; k <= 3; k++) {
                i2 = hb_linkage[labs(allobj[2][j])][k];
                fprintf(fp, format, xyz[i2][1], xyz[i2][2]);
            }
            fprintf(fp, "\n  gsave Ds LN grestore\n");
        }
    }
    close_file(fp);
    fp = open_file(temp_eps, "r");
    if (Gvars.label_RC8_YC6) {
        FILE *fp_label;
        fp_label = open_file(temp_labels, "w");
        while (fgets(str, sizeof str, fp) != NULL) {
            if (strstr(str, ") SCENTER"))
                fprintf(fp_label, "%s", str);
            else
                fprintf(fp_ok, "%s", str);
        }
        close_file(fp_label);
        fp_label = open_file(temp_labels, "r");
        while (fgets(str, sizeof str, fp_label) != NULL)
            fprintf(fp_ok, "%s", str);
        close_file(fp_label);
    } else {
        while (fgets(str, sizeof str, fp) != NULL)
            fprintf(fp_ok, "%s", str);
    }
    close_file(fp);
    fprintf(fp_ok, "\nshowpage\n");
}

static void stack2r3d(long *aidx, long **allobj, long num, double **xyz, long *ResSeq,
                      long nbond, long **linkage, long num_ring, long **ring_atom,
                      long num_hbond, long **hb_linkage, long **pair_num, long num_ball,
                      long num_residue, long *ibase, char *bseq, double scale_factor,
                      double ballrad, long *opts, FILE *fp)
{
    char label_style[BUF512];
    long i, ia, ib, j;
    long is_color, label_ring, frame_box;
    double dw, width3[4], hb_col[5], **atom_col, **base_col, **colrgb;
    label_ring = opts[3];
    is_color = opts[4];
    frame_box = opts[5];
    atom_col = dmatrix(0, NATOMCOL, 1, 3);
    base_col = dmatrix(0, NBASECOL, 1, 3);
    raster3d_header(num, xyz, scale_factor, 0, frame_box, fp);
    get_r3dpars(base_col, hb_col, width3, atom_col, label_style);
    if (!is_color)
        for (i = 1; i <= num_residue; i++)
            ibase[i] = NON_WC_IDX;
    if (num_ball) {
        colrgb = dmatrix(1, num_ball, 1, 3);
        for (i = 1; i <= num_ball; i++)
            (is_color) ? cpxyz(atom_col[aidx[i]], colrgb[i]) :
                cpxyz(base_col[NON_WC_IDX], colrgb[i]);
        fprintf(fp, "###\n### The following section is for atom sphere\n");
        cpk_model(num_ball, aidx, xyz, ballrad, colrgb, fp);
        free_dmatrix(colrgb, 1, num_ball, 1, 3);
    }
    fprintf(fp, "###\n### The following section is for bond-linkage\n");
    for (i = 1; i <= nbond; i++) {
        ia = linkage[i][1];
        ib = linkage[i][2];
        j = allobj[2][i + num_ball];
        dw = (lval_in_set(j, 1, pair_num[2][0], pair_num[2])) ? width3[3] : width3[2];
        r3d_rod(3, xyz[ia], xyz[ib], dw, base_col[ibase[j]], fp);
    }
    if (num_ring) {
        fprintf(fp, "###\n### The following section is for filled base-ring\n");
        fill_base_ring(num_residue, num_ring, ring_atom, xyz, ibase, bseq, base_col,
                       label_style, label_ring, ResSeq, fp);
    }
    if (num_hbond) {
        fprintf(fp, "###\n### The following section is for H-bonds\n");
        if (!is_color)
            cpxyz(base_col[NON_WC_IDX], hb_col);
        for (i = 1; i <= num_hbond; i++) {
            dw = (hb_linkage[i][1] == 2) ? width3[3] : width3[2];
            r3d_dash(xyz[hb_linkage[i][2]], xyz[hb_linkage[i][3]], dw, hb_col, fp);
        }
    }
    free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
    free_dmatrix(base_col, 0, NBASECOL, 1, 3);
}

static void get_stack_objs(long num_ball, double **xyz, long **seidx, long nbond,
                           long **linkage, long num_residue, long **ring_atom,
                           long num_hbond, long **hb_linkage, long **allobj)
{
    long i, i1, i2, iobj = 0, j, p1, p2;
    double dsum;
    for (i = 1; i <= num_ball; i++) {
        iobj++;
        allobj[1][iobj] = 0;
        allobj[2][iobj] = i;
        allobj[3][iobj] = lround(1000.0 * xyz[i][3]);
    }
    for (i = 1; i <= nbond; i++) {
        iobj++;
        allobj[1][iobj] = i;
        p1 = linkage[i][1];
        i1 = 0;
        for (j = 1; j <= num_residue; j++)
            if (lval_in_range(p1, seidx[j][1], seidx[j][2])) {
                i1 = j;
                break;
            }
        p2 = linkage[i][2];
        i2 = 0;
        for (j = 1; j <= num_residue; j++)
            if (lval_in_range(p2, seidx[j][1], seidx[j][2])) {
                i2 = j;
                break;
            }
        if (i1 != i2)
            allobj[2][iobj] = (i1 < i2) ? i1 : i2;
        else
            allobj[2][iobj] = i1;
        allobj[3][iobj] = lround(500.0 * (xyz[p1][3] + xyz[p2][3]));
    }
    for (i = 1; i <= num_residue; i++) {
        if (ring_atom[i][10] <= 0)
            continue;
        iobj++;
        allobj[1][iobj] = -1;
        allobj[2][iobj] = i;
        dsum = 0.0;
        for (j = 1; j <= ring_atom[i][10]; j++)
            dsum += xyz[ring_atom[i][j]][3];
        allobj[3][iobj] = lround(1000.0 * dsum / ring_atom[i][10]);
    }
    for (i = 1; i <= num_hbond; i++) {
        iobj++;
        allobj[1][iobj] = -2;
        allobj[2][iobj] = (hb_linkage[i][1] == 1) ? i : -i;
        allobj[3][iobj] = lround(500.0 * (xyz[hb_linkage[i][2]][3] + xyz[hb_linkage[i][3]][3]));
    }
}

static void process_stack(struct_args *args, long *opts)
{
    char *bseq, *ChainID, **AtomName, **ResName, **Miscs;
    double temp, **xyz;
    long i, nbond = 0, nbond_estimated, nobj, num, num_ball, num_residue;
    long bond_lmt, num_hbond = 0, num_ring = 0, **pair_num;
    long default_size[2] = { PS_DFTSIZE, FIG_DFTSIZE };
    long urxy[3], *RY, *aidx, *depth, *idx, *ibase, *ResSeq;
    long **allobj, **hb_linkage, **linkage, **ring_atom, **seidx;
    FILE *fp;
    num = number_of_atoms(args->inpfile, 1, "*");
    if (!num)
        fatal("PDB file without ATOM/HETATM records\n");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(args->inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    aidx = lvector(1, num);
    atom_idx(num, AtomName, Miscs, aidx);
    pair_num = lmatrix(1, 2, 0, MBASES);
    get_pair_num(args->inpfile, pair_num);
    if (opts[6])
        for (i = 1; i <= num; i++) {
            temp = xyz[i][1];
            xyz[i][1] = -xyz[i][2];
            xyz[i][2] = temp;
        }
    if (opts[7])
        get_side_view(1, num, xyz);
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);
    ring_atom = lmatrix(1, num_residue, 1, 10);
    if (opts[2])
        all_bring_atoms(num_residue, RY, seidx, AtomName, &num_ring, ring_atom);
    hb_linkage = lmatrix(1, num, 1, 3);
    if (opts[1] || opts[8]) {
        if (!strcmp(args->hbfile, "")) {
            if (pair_num[1][0] || pair_num[2][0])
                hbond_info(pair_num, bseq, seidx, aidx, AtomName, ResName, ChainID,
                           ResSeq, Miscs, xyz, RY, &num_hbond, hb_linkage);
            else
                hbond_pdb(num, num_residue, bseq, seidx, aidx, RY, AtomName, ResName,
                          ChainID, ResSeq, Miscs, xyz, &num_hbond, hb_linkage, opts[8]);
        } else
            read_hbinfo(args->hbfile, num, &num_hbond, hb_linkage);
    }
    bond_lmt = lround(NBOND_FNUM * num);
    nbond_estimated = 5 * bond_lmt;
    linkage = lmatrix(1, nbond_estimated, 1, 2);
    if (!strcmp(args->lkgfile, "")) {
        atom_linkage(1, num, aidx, xyz, Miscs, ChainID, nbond_estimated, &nbond, linkage);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, SNUM_FILE);
        write_lkglist(nbond, linkage, AtomName, ResName, ChainID, ResSeq, Miscs);
    } else
        read_lkginfo(args->lkgfile, num, &nbond, linkage);
    if (nbond > bond_lmt)
        fprintf(stderr, "too many linkages: *** possibly overlap ***\n");
    write_alc(num, nbond, AtomName, xyz, NULL, linkage, ATOMALC_FILE);
    ibase = lvector(1, num_residue);
    base_idx(num_residue, bseq, ibase, 0);
    num_ball = (args->ballrad == 0) ? 0 : num;
    nobj = num_ball + nbond + num_ring + num_hbond;
    allobj = lmatrix(1, 3, 1, nobj);
    get_stack_objs(num_ball, xyz, seidx, nbond, linkage, num_residue, ring_atom,
                   num_hbond, hb_linkage, allobj);
    fp = open_file(args->outfile, "w");
    if (opts[0] == 2)
        stack2r3d(aidx, allobj, num, xyz, ResSeq, nbond, linkage, num_ring,
                  ring_atom, num_hbond, hb_linkage, pair_num, num_ball,
                  num_residue, ibase, bseq, args->scale_factor, args->ballrad, opts, fp);
    else {
        idx = lvector(1, nobj);
        lsort(nobj, allobj[3], idx);
        adjust_xy(num, xyz, 0L, NULL, args->scale_factor, default_size[opts[0]], urxy);
        depth = lvector(1, nobj);
        if (opts[0] == 1) {
            get_depth(nobj, allobj[3], depth);
            get_fig_xy(num, xyz, 0L, NULL, urxy, opts[5], fp);
            stack2fig(nobj, idx, aidx, depth, allobj, xyz, linkage, hb_linkage,
                      ring_atom, ibase, opts[4], pair_num, opts[3], ResSeq, args->ballrad, fp);
        } else {
            get_ps_xy(args->outfile, urxy, opts[5], fp);
            stack2ps(nobj, idx, aidx, allobj, xyz, linkage, hb_linkage, ring_atom,
                     ibase, opts[4], pair_num, opts[3], ResSeq, args->ballrad, fp);
        }
        free_lvector(idx, 1, nobj);
        free_lvector(depth, 1, nobj);
    }
    close_file(fp);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(pair_num, 1, 2, 0, MBASES);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_lmatrix(ring_atom, 1, num_residue, 1, 10);
    free_lmatrix(hb_linkage, 1, num, 1, 3);
    free_lvector(ibase, 1, num_residue);
    free_lmatrix(linkage, 1, nbond_estimated, 1, 2);
    free_lvector(aidx, 1, num);
    free_lmatrix(allobj, 1, 3, 1, nobj);
}

int main(int argc, char *argv[])
{
    struct_args args;
    long opts[BUF512] = { 0 };
    set_my_globals(argv[0]);
    stack2img_cmdline(argc, argv, &args, opts);
    if (strcmp(args.hbfile, HB_FILE))
        remove_file(HB_FILE);
    if (strcmp(args.lkgfile, LKG_FILE))
        remove_file(LKG_FILE);
    remove_file(ATOMALC_FILE);
    remove_file(SNUM_FILE);
    if (Gvars.label_RC8_YC6)
        opts[3] = TRUE;
    process_stack(&args, opts);
    clear_my_globals();
    return 0;
}
