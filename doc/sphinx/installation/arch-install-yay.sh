#!/bin/bash

mkdir -p ~/yay
cd ~/yay
git clone https://aur.archlinux.org/yay.git .
makepkg -si
