module load singularity/3.8.4

export SIMG_IMAGE=/packages/containers/Rstudio/rstudio-with_modules-4.2.0-2.sif

export RSLIB=$HOME/.rstudio-server/lib
export RSRUN=$HOME/.rstudio-server/run
export RSTMP=/scratch/${USER}/.rstudio-server/tmp

singularity shell -B /labs/banovich -B /scratch/${USER} -B /packages/petagene/protect_1.3.21-1/bin -B /opt,$RSTMP:/tmp,$RSLIB:/var/lib/rstudio-server/,$RSRUN:/var/run/rstudio-server/ $SIMG_IMAGE