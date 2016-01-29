# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
export HISTCONTROL=ignoreboth:erasedups  # no duplicate entries
shopt -s histappend
export HISTSIZE=9999
export HISTFILESIZE=999999

export PATH=$PATH:/share/PI/nburns/mawk-1.3.4-20150503:/share/PI/nburns/orca:/share/PI/nburns/cfour_v1/bin

ulimit -c 0
module load gaussian/g09
module load openmpi/1.8.7

alias squeue1="squeue -u ekwan16 -o \"%.10i %.15P  %.50j %.15u %.8T %.10M %.15l %.6D %.6m %.2C %R\""
alias squeue2="squeue -o \"%.10i %.15P  %.50j %.15u %.8T %.10M %.15l %.6D %.6m %.2C %R\""
alias ls="ls --color=auto --group-directories-first"

export PS1='[\[\033[01;32m\]\u\[\033[00m\]@\h \w]\$ '
