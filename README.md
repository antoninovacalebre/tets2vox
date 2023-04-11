# tets2vox
Converts tetrahedral meshes (MSH 2.2) and voxel meshes (VoxHenry) to voxel meshes (MSH 2.2, VoxHenry).

## Usage

Input file format: MSH 2.2 tetrahedra (`.msh` or `.msh2`), VoxHenry voxels (`.vhr`)
Output file format: MSH 2.2 voxels (`.msh` or `.msh2`), VoxHenry voxels (`.vhr`)

```
Usage: tets2vox.exe [OPTIONS] --input <input> --output <output>

Options:
  -i, --input <input>    Input file
  -o, --output <output>  Output file
  -r, --res <res>        Resolution (default 100)
  -h, --help             Print help
  -V, --version          Print version
```
