FROM jupyter/scipy-notebook
MAINTAINER Dave Williams <cdave@uw.edu>

USER root

## Geometric operation dependencies 
RUN apt-get update && apt-get install -y \
  blender \
  openscad 

## Configure Jupiter notebook extensions for nice environment
    # Enable JS widgets
    RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension
    # Install contributed notebook extensions
    RUN pip install https://github.com/ipython-contrib/jupyter_contrib_nbextensions/tarball/master
    RUN jupyter contrib nbextension install --system
    RUN jupyter nbextension enable collapsible_headings/main
    RUN jupyter nbextension enable spellchecker/main
    # Install vim bindings
    RUN jupyter nbextension install https://github.com/lambdalisue/jupyter-vim-binding/archive/master.tar.gz --system
    RUN jupyter nbextension enable jupyter-vim-binding-master/vim_binding


# Install python dependencies and package
WORKDIR $HOME/celloutline
COPY . .
RUN sudo chown -R $NB_USER:$NB_GID . 
USER $NB_USER
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install -e .
WORKDIR $HOME

