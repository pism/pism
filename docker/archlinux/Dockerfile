FROM archlinux/base:latest
RUN pacman -Sy && pacman --noconfirm -S git base-devel

RUN useradd --create-home --shell=/bin/false builder && usermod --lock builder
RUN echo "builder ALL = NOPASSWD: /usr/bin/pacman" >> /etc/sudoers
USER builder

WORKDIR /home/builder/
RUN git clone https://aur.archlinux.org/yay.git && cd yay && makepkg --noconfirm -si

RUN yay --noconfirm --answerdiff=None -Sy pism
