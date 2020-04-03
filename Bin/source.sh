
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/nfs1/public2/User/develop/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/nfs1/public2/User/develop/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/nfs1/public2/User/develop/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/nfs1/public2/User/develop/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate TheThinker
# <<< conda initialize <<<

