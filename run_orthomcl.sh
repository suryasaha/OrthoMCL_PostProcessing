#!/bin/sh

# Surya Saha
# Purpose: Run orthomcl commands for each parameter setting

set -o nounset
set -o errexit

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))

usage() {
  cat << EOF
    usage:
    $PROGNAME <orthomcl install directory> <config file for parameters>
    
    Example:
    $PROGNAME /home/surya/tools/orthomclSoftware-v2.0.9 orthomcl.config.eval-5percent30
       
EOF
}

printf "Make sure your protein FAA files are in Refseq format\nExample: >gi|16262454|ref|NP_435247.1| FdoG formate dehydrogenase-O alpha subunit [Sinorhizobium meliloti 1021]\nUse orthomcl.convert2RefseqFasta.pl for formatting.\n"


if [ "$#" -ne 2 ]
then
	usage
fi

printf "$1\n"
printf "$2\n"
#${1}/bin/orthomclInstallSchema

${1}/bin/orthomclInstallSchema "$2" ${2}.install_log
${1}/bin/orthomclLoadBlast "$2" similarSequences.txt
${1}/bin/orthomclPairs "$2" ${2}.run_log cleanup=yes
${1}/bin/orthomclDumpPairsFiles "$2" 

SUFFIX=`echo $2| sed 's,orthomcl.config.,,'`
#printf $SUFFIX
mv mclInput mclInput.${SUFFIX}
mv pairs/ pairs-${SUFFIX}

mcl mclInput.${SUFFIX} --abc -I 1.5 -o mclOutput.${SUFFIX}
${1}/bin/orthomclMclToGroups Clostr 1000 < mclOutput.${SUFFIX} > groups.txt.${SUFFIX}

