#include "samtools.pysam.h"

/*  dict.c -- create a sequence dictionary file.

    Copyright (C) 2015, 2020 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <getopt.h>
#include "htslib/kseq.h"
#include "htslib/hts.h"

KSEQ_INIT(gzFile, gzread)

typedef struct _args_t
{
    char *output_fname, *fname;
    char *assembly, *species, *uri;
    int  alias, header;
}
args_t;

static void write_dict(const char *fn, args_t *args)
{
    hts_md5_context *md5;
    int l, i, k;
    gzFile fp;
    kseq_t *seq;
    unsigned char digest[16];
    char hex[33];

    fp = strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        fprintf(samtools_stderr, "dict: %s: No such file or directory\n", fn);
        samtools_exit(1);
    }
    FILE *out = samtools_stdout;
    if (args->output_fname) {
        out = fopen(args->output_fname, "w");
        if (out == NULL) {
          fprintf(samtools_stderr, "dict: %s: Cannot open file for writing\n", args->output_fname);
          samtools_exit(1);
        }
    }

    if (!(md5 = hts_md5_init()))
        samtools_exit(1);

    seq = kseq_init(fp);
    if (args->header) fprintf(out, "@HD\tVN:1.0\tSO:unsorted\n");
    while ((l = kseq_read(seq)) >= 0) {
        for (i = k = 0; i < seq->seq.l; ++i) {
            if (seq->seq.s[i] >= '!' && seq->seq.s[i] <= '~')
                seq->seq.s[k++] = toupper(seq->seq.s[i]);
        }
        hts_md5_reset(md5);
        hts_md5_update(md5, (unsigned char*)seq->seq.s, k);
        hts_md5_final(digest, md5);
        hts_md5_hex(hex, digest);
        fprintf(out, "@SQ\tSN:%s\tLN:%d\tM5:%s", seq->name.s, k, hex);
        if (args->alias) {
            const char *name = seq->name.s;
            if (strncmp(name, "chr", 3) == 0) {
                name += 3;
                fprintf(out, "\tAN:%s", name);
            }
            else
                fprintf(out, "\tAN:chr%s", name);

            if (strcmp(name, "M") == 0)
                fprintf(out, ",chrMT,MT");
            else if (strcmp(name, "MT") == 0)
                fprintf(out, ",chrM,M");
        }
        if (args->uri)
            fprintf(out, "\tUR:%s", args->uri);
        else if (strcmp(fn, "-") != 0) {
#ifdef _WIN32
            char *real_path = _fullpath(NULL, fn, PATH_MAX);
#else
            char *real_path = realpath(fn, NULL);
#endif
            fprintf(out, "\tUR:file://%s", real_path);
            free(real_path);
        }
        if (args->assembly) fprintf(out, "\tAS:%s", args->assembly);
        if (args->species) fprintf(out, "\tSP:%s", args->species);
        fprintf(out, "\n");
    }
    kseq_destroy(seq);
    hts_md5_destroy(md5);

    if (args->output_fname) fclose(out);
    gzclose(fp);
}

static int dict_usage(void)
{
    fprintf(samtools_stderr, "\n");
    fprintf(samtools_stderr, "About:   Create a sequence dictionary file from a fasta file\n");
    fprintf(samtools_stderr, "Usage:   samtools dict [options] <file.fa|file.fa.gz>\n\n");
    fprintf(samtools_stderr, "Options: -a, --assembly STR    assembly\n");
    fprintf(samtools_stderr, "         -A, --alias, --alternative-name\n");
    fprintf(samtools_stderr, "                               add AN tag by adding/removing 'chr'\n");
    fprintf(samtools_stderr, "         -H, --no-header       do not print @HD line\n");
    fprintf(samtools_stderr, "         -o, --output FILE     file to write out dict file [samtools_stdout]\n");
    fprintf(samtools_stderr, "         -s, --species STR     species\n");
    fprintf(samtools_stderr, "         -u, --uri STR         URI [file:///abs/path/to/file.fa]\n");
    fprintf(samtools_stderr, "\n");
    return 1;
}

int dict_main(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->header = 1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"no-header", no_argument, NULL, 'H'},
        {"alias", no_argument, NULL, 'A'},
        {"alternative-name", no_argument, NULL, 'A'},
        {"assembly", required_argument, NULL, 'a'},
        {"species", required_argument, NULL, 's'},
        {"uri", required_argument, NULL, 'u'},
        {"output", required_argument, NULL, 'o'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ( (c=getopt_long(argc,argv,"?AhHa:s:u:o:",loptions,NULL))>0 )
    {
        switch (c)
        {
            case 'A': args->alias = 1; break;
            case 'a': args->assembly = optarg; break;
            case 's': args->species = optarg; break;
            case 'u': args->uri = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'H': args->header = 0; break;
            case 'h': return dict_usage();
            default: return dict_usage();
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(STDIN_FILENO) ) fname = "-";  // reading from stdin
        else return dict_usage();
    }
    else fname = argv[optind];

    write_dict(fname, args);
    free(args);
    return 0;
}
