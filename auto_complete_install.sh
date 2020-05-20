#!/usr/bin/env bash
if [[ $(basename $SHELL) = 'bash' ]];
then
    if [ -f ~/.bashrc ];
    then
        echo "Installing bash autocompletion..."
        grep -q '_NOMAD_COMPLETE' ~/.bashrc
        if [[ $? -ne 0 ]]; then
            echo "" >> ~/.bashrc
            echo 'eval "$(_NOMAD_COMPLETE=source nomad)"' >> ~/.bashrc
        fi
    fi
elif [[ $(basename $SHELL) = 'zsh' ]];
then
    if [ -f ~/.zshrc ];
    then
        echo "Installing zsh autocompletion..."
        grep -q 'nomad-autocompletion' ~/.zshrc
        if [[ $? -ne 0 ]]; then
            echo "" >> ~/.zshrc
            echo "autoload bashcompinit" >> ~/.zshrc
            echo "bashcompinit" >> ~/.zshrc
            echo 'eval "$(_NOMAD_COMPLETE=source nomad)"' >> ~/.nomad-autocompletion.sh
            echo "source ~/.nomad-autocompletion.sh" >> ~/.zshrc
        fi
    fi
fi