# Bionano_LocalAssembler
Bionano local comparative assembler for complex regions. Experimental version

usage: cluster.py [-h] [-x XMAP] [-q QMAP] [-r RMAP] [-o OUT] [-c CHR]

                  [-s START] [-e END] [-u UPSTREAM_EXTEND]

                  [-d DOWNSTREAM_EXTEND] [-m CLUSTER_MAX]

                  [-p STDSHIFT_PENALTY] [-b BASE_PENALTY]



script to cluster and assemble bionano local alignments


optional arguments:

  -h, --help            show this help message and exit

  -x XMAP, --xmap XMAP  xmap file

  -q QMAP, --query_cmap QMAP

                        query_cmap file

  -r RMAP, --ref_cmap RMAP

                        rmap file

  -o OUT, --output OUT  output file

  -c CHR, --chr CHR     target region chromosome

  -s START, --start START

                        target region start

  -e END, --end END     target region end

  -u UPSTREAM_EXTEND, --upstream_extend UPSTREAM_EXTEND

                        extending anchor region upstream size

  -d DOWNSTREAM_EXTEND, --downstream_extend DOWNSTREAM_EXTEND

                        extending anchor region downstream size

  -m CLUSTER_MAX, --cluster_max CLUSTER_MAX

                        maximum difference to cluster

  -p STDSHIFT_PENALTY, --stdshift_penalty STDSHIFT_PENALTY

                        strandshift_penalty

  -b BASE_PENALTY, --base_penalty BASE_PENALTY

                        defaulty baseline penalty





It finds all alignments within anchor(target) region plus the extend region. it will also extend corresponding size on the query molecules. eg. if the most right alignment is 14647178+700000 on the reference, it will extend 750000-700000 downstream. 

-d is the donestream size, and
-u is the upstream size. 

their defaults are zeros, which means no extending. 
