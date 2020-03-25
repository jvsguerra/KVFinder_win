# KVFinder_win
Deprecated KVFinder version

Working features:
- Command-line interface (kvfinder-win.exe)
- Graphical user interface (kvfinder_pymolplugin_windows.py)

Compatible with:
- Windows 10
- PyMOL 1.8x (check *Installation of PyMOL 1.8x* section)
- Python 2

# Download and Installation
- Get the latest release from [here](https://github.com/jvsguerra/KVFinder_win/releases/).
- Unzip the KVFinder_win-x.x.zip.
- In the environment variables, add KVFinder_PATH variable as KVFinder-win-x.x directory location.

# Install KVFinder PyMOL plugin
- Open PyMOL 1.8x
- Go to the **Plugin Manager** option under **Plugin** menu
- Go to the **Install New Plugin**, select **Choose file...** inside **Install from local file**.
- Select **kvfinder_pymolplugin_windows.py** inside **tools** directory.
- A window will appear confirming that the plugin has been installed.
- Restart PyMOL.
- KVFinder plugin is ready to use.

# Installation of PyMOL 1.8x
- Follow steps from http://tubiana.me/how-to-install-and-compile-pymol-windows-linux-mac/ 

# Example of use
We include a example directory with a PDB (1FMO.pdb) without water and its ligands.

To run KVFinder with command-line interface (kvfinder-win.exe):
```
# Inside KVFinder-win directory
C:\Users\User\path_to_KVFinder\KVFinder-win-x.x> kvfinder-win.exe -h
KVFinder Help
-i PDB File for a single run
-l File with a list of PDBs for several runs
-p Parameters File
-g Define grid size
-r Probe In size
-t Probe Out size
-v Volume filter

C:\Users\User\path_to_KVFinder\KVFinder-win-x.x> kvfinder-win.exe -i example\1FMO.pdb
Running KVFinder for input\1FMO.pdb...done!

```

# Citing KVFinder software
If you use KVFinder in your publication, refer to this [article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-197).

Oliveira, S.H., Ferraz, F.A., Honorato, R.V. et al. KVFinder: steered identification of protein cavities as a PyMOL plugin. BMC Bioinformatics 15, 197 (2014). https://doi.org/10.1186/1471-2105-15-197

