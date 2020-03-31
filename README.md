# CTC model
This repository includes the necessary data, scripts and python code to reproduce the results of:

F Font-Clos, S Zapperi, CAM La Porta. *Humanoid global vascular model for circulating tumor cells*, The Journal, 1234-56, 2020.  
[link to publication](http://add-the-link-when-accepted)

In addition, we provide code to solve hemodynamic equations in a directed network, simulate Circulating Tumor Cells (CTC) trajectories on it, and render a 3D-like movie of CTC trajectories in pancreatic cancer. You can see the final video [here](http://add-link-to-youtube-or-vimeo). 

## Simulating CTC trajectories
You can simulate your own CTC trajectories using the `code/launch_tracers_fullsystem.py` file:
```
python code/launch_tracers_fullsystem.py --help
usage: Launch many tracers from an organ into the blood stream
       --organ ORGAN --network NETWORK -n NUM_TRACERS --output_dir OUTPUT_DIR
       [--verbose] [--drop_trajectories] [--keep_times] [--suffix SUFFIX] [-h]

required arguments:
  --organ ORGAN         Name of recognized organ. See Makefile for list of
                        available organs
  --network NETWORK     Full body network with attachment probabilities
  -n NUM_TRACERS, --num_tracers NUM_TRACERS
                        Number of tracers (to be splitted among available
                        starting points
  --output_dir OUTPUT_DIR
                        Path to output dir

optional arguments:
  --verbose             This is a binary switch flag.
  --drop_trajectories   do not store full trajectories, just final position
  --keep_times          additionaly store simulation time
  --suffix SUFFIX       Suffix to add to output files
  -h, --help            Show this help message and exit.
```
For instance, to simulate 100 CTC trajectories starting from the pancreas, and store the output in certain `path/to/output/dir`, labeling your simulations as `_my_trajectories`, and keeping all information (will produce larger output), you would run:
```bash
python code/launch_tracers_fullsystem.py \
				--keep_times \
				--organ pancreas \
				--map_endpoints \
				--suffix _my_trajectories \
				--network output/data/solved_flow_networks/solved_total_closed.p \
				--num_tracers 100 \
				--output_dir path/to/output/dir
```
This will produce the file `path/to/output/dir/pancreas_tracers_100_my_trajectories.p` with a size of around 300Mb. The output is stored as a pickled pandas DataFrame, and includes all positions, times and labels of the associated network nodes, for all steps of 100 CTC simulated trajectories. Notice that if you want to compute organ-to-organ probabilities, you can pass the `--drop_trajectories` flag and omit the `--keep_times` one to produce much smaller output files. 

## Reproducing the figures of the manuscript
You can exactly reproduce the manuscript figures by executing the notebooks in the `notebooks/` directory:
  1. [Solve hemodynamic flow eqs on network](notebooks/1-solve-hemodynamic-flow-eqs-on-network.ipynb)
  2. [Compute fraction of cells from simulations](notebooks/2-compute-fraction-of-cells-from-simulations.ipynb)
  3. [Autopsy metastasis freqs](notebooks/3-autopsy-metastasis-freqs.ipynb)
  4. [Model fraction cells heatmap](notebooks/4-model-fraction-cells-heatmap.ipynb)
  5. [Whole-body flow plots](notebooks/5-whole-body-flow-plots.ipynb)
  6. [Model vs data plots](notebooks/6-model-data-plots.ipynb)
  7. [Trajectories 2D pancreas](notebooks/7-trajectories-2d-pancreas-liver.ipynb)

Alternatively, you can run all notebooks in one go with 
```bash
./run_all_notebooks.sh
```
which will re-generate all figures of the manuscript.

## Rendering a 3D simulation of pancreatic cancer CTC
To render the CTC movie from scratch, you need to install `blender 2.8` on a high-performance computer. Notice the video is rendered at very high resolution. The scripts provided in `video/` can be used generate 1250 PNG images of size 3840 x 2160 pixels, each of which has a size of 5-15Mb. We include a few sample frames in the `video/frames` directory.

### Installing blender
To install blender on Ubuntu machines, you can use the `video/install_blender.sh` script. Notice you will need administrator privileges for that. Simply run
```bash
cd video
./install_blender.sh
```

### Getting the blender file
Our blender raw file `video/ctc_movie.blend` is too large to be hosted on github, but you can get it via this google drive link: [ctc_movie.blend](https://drive.google.com/open?id=1b95N8FhRqwPo0G7Q-YjtbZiHUNf1w1cE). Before trying to render frames, download the blender file and place it in the `video/` folder.

### Rendering all frames
**Note:** we recommend rendering on a large cluster, AWS instance, renderfarm service or similar. Rendering on a standard desktop computer is possible but might weeks.

To render all frames, run the `video/install_blender.sh` script. 
```bash
cd video
./render_frames.sh
```

### Encoding the video
We use `ffmpeg` to encode the video from individual PNG frames. Once you have rendered all frames, to obtain large, medium and small size videos in one go run
```bash
cd video
./make_ctc_movie.sh
```
## Data sources
The original **mesh files** from which the flow network is obtained have been downloaded from the BodyParts3D project:
+ BodyParts3D/Anatomography: [link to main site](https://lifesciencedb.jp/bp3d/) | [link to download site](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html)

The **autopsy data** comes from several publications, see the manuscript for details. We include hand-curated excel files from which it is easy to retrieve the data used in the manuscript:
+ Abrams, 1950: [publication](http://dx.doi.org/10.1002/1097-0142(1950)3:1%3C74::aid-cncr2820030111%3E3.0.co;2-7) | [excel file](data/metastasis_frequency_data/Abrams1950/Abrams1950-longform-byhand.xlsx)
+ Bubendorf, 2000: [publication](http://dx.doi.org/10.1053/hp.2000.6698) | [excel file](data/metastasis_frequency_data/Bubendorf2000/Bubendorf2000-longform-byhand.xlsx)
+ Budczies, 2015: [publication](http://dx.doi.org/10.18632/oncotarget.2677) | [excel file](data/metastasis_frequency_data/Budczies2015/Budczies2015-longform-byhand.xlsx)
+ diSibio, 2008: [publication](http://dx.doi.org/10.1043/1543-2165(2008)132[931:MPOCRF]2.0.CO;2) | [excel file](data/metastasis_frequency_data/diSibio2008/diSibio2008-longform-byhand.xlsx)
+ Schlageter, 2016: [publication](http://dx.doi.org/10.1159/000446245) | [excel file](data/metastasis_frequency_data/Schlageter2016/Schlageter2016-longform-byhand.xlsx)


## Dependencies
Beyond well-known python packages (numpy, matplotlib, seaborn, pandas, etc), this repository depends on the following: 
+ [pyembree](https://github.com/scopatz/pyembree)
+ [trimesh](https://trimsh.org/) (forked version [here](https://github.com/fontclos/trimesh))
+ [K3D-jupyter](https://github.com/K3D-tools/K3D-jupyter)
+ [blender (2.8)](https://www.blender.org/)(installation recommended via this [install script for ubuntu](video/install_blender.sh))

## Setup (conda)
We recommend installing dependencies in a clean, isolated environment via `pipenv`, `virtualenv` or similar. For instance, if you are a `conda` user, you could do:
```bash
conda create --yes --name ctcmodel pip
conda activate ctcmodel
pip install -r requirements.txt
```

