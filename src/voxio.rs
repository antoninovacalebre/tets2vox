use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use ndarray::prelude::*;

pub fn vhr2vox(file: &str) -> (Array3<i64>, f64, Array1<f64>) {
    // VHR file format :
    // - Commented lines start with % or *
    // - we can ignore the ports

    // Example :
    // freq=1.0
    // dx=0.0001
    // LMN=100,100,100
    //
    // StartVoxelList
    // V 0 0 0 1.000000e+00
    // V 0 0 1 1.000000e+00
    // ...
    // EndVoxelList
    // N port1 N 0 0 0 -z
    // N port2 N 99 99 99 z

    let mut vhrvoxels: HashMap<(i64, i64, i64), i64> = HashMap::new();
    let mut dx: f64 = 0.0;
    let mut lmn: [usize; 3] = [0, 0, 0];

    #[derive(Debug)]
    enum VHRState {
        Preamble,
        VoxelList,
        PortList,
    }

    let mut state = VHRState::Preamble;

    let lines = std::fs::read_to_string(file).unwrap();
    let lines: Vec<&str> = lines.lines().collect();

    for line in lines {
        let line = line.trim();
        
        if line.starts_with("%") || line.starts_with("*") {
            continue;
        }
        
        match state {
            VHRState::Preamble => {
                if line.starts_with("dx=") {
                    dx = line.split("=").collect::<Vec<&str>>()[1]
                        .trim()
                        .parse::<f64>()
                        .unwrap();
                } else if line.starts_with("LMN=") {
                    let lmn_str = line.split("=").collect::<Vec<&str>>()[1].trim();
                    let lmn_str = lmn_str.split(",").collect::<Vec<&str>>();
                    lmn[0] = lmn_str[0].parse::<usize>().unwrap();
                    lmn[1] = lmn_str[1].parse::<usize>().unwrap();
                    lmn[2] = lmn_str[2].parse::<usize>().unwrap();
                } else if line.starts_with("StartVoxelList") {
                    state = VHRState::VoxelList;
                } else {
                    continue;
                }
            }
            VHRState::VoxelList => {
                if line.starts_with("EndVoxelList") {
                    state = VHRState::PortList;
                } else if line.starts_with("V") {
                    let line = line.split_whitespace().collect::<Vec<&str>>();
                    let i = line[1].parse::<i64>().unwrap() - 1;
                    let j = line[2].parse::<i64>().unwrap() - 1;
                    let k = line[3].parse::<i64>().unwrap() - 1;
                    vhrvoxels.insert((i, j, k), 1);
                } else {
                    continue;
                }
            }
            VHRState::PortList => {
                continue;
            }
        };
    }

    let mut voxels = Array3::<i64>::zeros((lmn[0], lmn[1], lmn[2]));
    for ((i, j, k), _v) in vhrvoxels {
        voxels[[i as usize, j as usize, k as usize]] = 1;
    }

    (voxels, dx, array![0.0, 0.0, 0.0])
}

pub fn vox2vhr(voxels: &Array3<i64>, dx: f64, file: &str) {
    let mut gnd_vox: Array1<i64> = array![0, 0, 0];
    let mut gnd_vox_set: bool = false;

    let f = File::create(file).unwrap();
    let mut writer = BufWriter::new(f);

    writeln!(writer, "freq=1.0").unwrap();
    writeln!(writer, "dx={}", dx).unwrap();
    writeln!(
        writer,
        "LMN={},{},{}",
        voxels.shape()[0],
        voxels.shape()[1],
        voxels.shape()[2]
    )
    .unwrap();
    writeln!(writer, "").unwrap();

    writeln!(writer, "StartVoxelList").unwrap();
    for i in 0..voxels.shape()[0] {
        for j in 0..voxels.shape()[1] {
            for k in 0..voxels.shape()[2] {
                if voxels[[i, j, k]] > 0 {
                    if !gnd_vox_set {
                        gnd_vox = array![i as i64, j as i64, k as i64];
                        gnd_vox_set = true;
                    }

                    writeln!(
                        writer,
                        "V {} {} {} {:e}",
                        i + 1,
                        j + 1,
                        k + 1,
                        voxels[[i, j, k]] as f64
                    )
                    .unwrap();
                }
            }
        }
    }
    writeln!(writer, "EndVoxelList").unwrap();
    writeln!(
        writer,
        "N port1 N {} {} {} -z",
        gnd_vox[0] + 1,
        gnd_vox[1] + 1,
        gnd_vox[2] + 1
    )
    .unwrap();

    writer.flush().unwrap();
}

pub fn vox2gmsh(voxels: &Array3<i64>, dx: f64, mesh_center: &Array1<f64>, file: &str) {
    // Map of faces to stitch together as a list
    // of tuples (vector_offsets, faces_to_stitch)
    let stitch_map: Vec<(Vec<i64>, Vec<i64>)> = vec![
        (vec![0, 0, 1], vec![-1, -1, 1, 0, -1, -1, 5, 4]), // top
        (vec![0, 0, -1], vec![3, 2, -1, -1, 7, 6, -1, -1]), // bot
        (vec![0, -1, 0], vec![-1, -1, -1, -1, 0, 1, 2, 3]), // front
        (vec![0, 1, 0], vec![4, 5, 6, 7, -1, -1, -1, -1]), // back
        (vec![1, 0, 0], vec![-1, 0, 3, -1, -1, 4, 7, -1]), // right
        (vec![-1, 0, 0], vec![1, -1, -1, 2, 5, -1, -1, 6]), // left
        (vec![-1, 1, -1], vec![6, -1, -1, -1, -1, -1, -1, -1]), // 0
        (vec![1, 1, -1], vec![-1, 7, -1, -1, -1, -1, -1, -1]), // 1
        (vec![1, 1, 1], vec![-1, -1, 4, -1, -1, -1, -1, -1]), // 2
        (vec![-1, 1, 1], vec![-1, -1, -1, 5, -1, -1, -1, -1]), // 3
        (vec![-1, -1, -1], vec![-1, -1, -1, -1, 2, -1, -1, -1]), // 4
        (vec![1, -1, -1], vec![-1, -1, -1, -1, -1, 3, -1, -1]), // 5
        (vec![1, -1, 1], vec![-1, -1, -1, -1, -1, -1, 0, -1]), // 6
        (vec![-1, -1, 1], vec![-1, -1, -1, -1, -1, -1, -1, 1]), // 7
        (vec![0, 1, -1], vec![7, 6, -1, -1, -1, -1, -1, -1]), // 0-1
        (vec![1, 1, 0], vec![-1, 4, 7, -1, -1, -1, -1, -1]), // 1-2
        (vec![0, 1, 1], vec![-1, -1, 4, 5, -1, -1, -1, -1]), // 2-3
        (vec![-1, 1, 0], vec![5, -1, -1, 6, -1, -1, -1, -1]), // 3-0
        (vec![0, -1, -1], vec![-1, -1, -1, -1, 3, 2, -1, -1]), // 4-5
        (vec![1, -1, 0], vec![-1, -1, -1, -1, -1, 0, 3, -1]), // 5-6
        (vec![0, -1, 1], vec![-1, -1, -1, -1, -1, -1, 1, 0]), // 6-7
        (vec![-1, -1, 0], vec![-1, -1, -1, -1, 1, -1, -1, 2]), // 7-4
        (vec![-1, 0, -1], vec![2, -1, -1, -1, 6, -1, -1, -1]), // 0-4
        (vec![1, 0, -1], vec![-1, 3, -1, -1, -1, 7, -1, -1]), // 1-5
        (vec![1, 0, 1], vec![-1, -1, 0, -1, -1, -1, 4, -1]), // 2-6
        (vec![-1, 0, 1], vec![-1, -1, -1, 1, -1, -1, -1, 5]), // 3-7
    ];

    let offsets = vec![
        vec![-0.5, 0.5, -0.5],
        vec![0.5, 0.5, -0.5],
        vec![0.5, 0.5, 0.5],
        vec![-0.5, 0.5, 0.5],
        vec![-0.5, -0.5, -0.5],
        vec![0.5, -0.5, -0.5],
        vec![0.5, -0.5, 0.5],
        vec![-0.5, -0.5, 0.5],
    ];
    // multiply all offsets by dx
    let offsets = offsets
        .iter()
        .map(|x| x.iter().map(|y| y * dx).collect())
        .collect::<Vec<Vec<f64>>>();

    let grid_center = vec![
        0.5 * (voxels.shape()[0] as f64) * dx,
        0.5 * (voxels.shape()[1] as f64) * dx,
        0.5 * (voxels.shape()[2] as f64) * dx,
    ];

    let ncell = voxels.iter().filter(|&&x| x != 0).count();

    let mut all_faces: Vec<Vec<i64>> = Vec::with_capacity(ncell * 6);
    let mut nodes = vec![];
    let mut cubes: HashMap<(usize, usize, usize), Vec<i64>> = HashMap::new();

    for (voxel, v) in voxels.indexed_iter() {
        if *v == 0 {
            continue;
        }

        // move mesh to origin as minimum and find cube
        let mut coord = vec![voxel.0 as f64, voxel.1 as f64, voxel.2 as f64];
        coord = coord.iter().map(|x| (x + 0.5) * dx).collect();

        // Initialize cube as 1D vector of 8 nodes (-1)
        let mut cube = Array::from_elem((8,), -1);

        // Iterate over stich_map
        for (_, (neighbour_offset, faces)) in stitch_map.iter().enumerate() {
            // check if tuple is a key in cubes
            let neighbour = (
                (voxel.0 as i64 + neighbour_offset[0]) as usize,
                (voxel.1 as i64 + neighbour_offset[1]) as usize,
                (voxel.2 as i64 + neighbour_offset[2]) as usize,
            );
            if cubes.contains_key(&neighbour) {
                // iterate over faces
                for (j, node) in faces.iter().enumerate() {
                    if *node >= 0 {
                        cube[j] = cubes[&neighbour][*node as usize];
                    }
                }
            }
        }

        // Iterate over cube
        for i in 0..cube.len() {
            let node = cube[i];
            if node < 0 {
                let mut new_node = vec![
                    coord[0] + offsets[i][0] - (grid_center[0] - mesh_center[0]),
                    coord[1] + offsets[i][1] - (grid_center[1] - mesh_center[1]),
                    coord[2] + offsets[i][2] - (grid_center[2] - mesh_center[2]),
                ];
                new_node = new_node.iter().map(|x| *x).collect();
                nodes.push(new_node);
                cube[i] = nodes.len() as i64 - 1;
            }
        }

        cubes.insert(voxel, cube.to_vec());

        all_faces.push(vec![cube[4], cube[5], cube[6], cube[7]]);
        all_faces.push(vec![cube[1], cube[0], cube[3], cube[2]]);
        all_faces.push(vec![cube[1], cube[5], cube[6], cube[2]]);
        all_faces.push(vec![cube[0], cube[4], cube[7], cube[3]]);
        all_faces.push(vec![cube[7], cube[6], cube[2], cube[3]]);
        all_faces.push(vec![cube[0], cube[1], cube[5], cube[4]]);
    }

    // Filter all faces that have no duplicates
    let bfaces = all_faces
        .into_iter() // convert to iterator
        .collect::<Vec<Vec<i64>>>()
        .into_iter()
        .collect::<std::collections::HashSet<Vec<i64>>>() // convert to hash set
        .into_iter()
        .collect::<Vec<Vec<i64>>>(); // convert back to vector

    // open file
    let f = File::create(file).unwrap();
    let mut writer = BufWriter::new(f);

    // write header
    writeln!(writer, "$MeshFormat").unwrap();
    writeln!(writer, "2.2 0 8").unwrap();
    writeln!(writer, "$EndMeshFormat").unwrap();

    // write nodes
    writeln!(writer, "$Nodes").unwrap();
    writeln!(writer, "{}", nodes.len()).unwrap();
    for (i, mnode) in nodes.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            i + 1,
            mnode[0],
            mnode[1],
            mnode[2]
        )
        .unwrap();
    }
    writeln!(writer, "$EndNodes").unwrap();

    // write elements
    writeln!(writer, "$Elements").unwrap();
    writeln!(writer, "{}", cubes.len() + bfaces.len()).unwrap();
    let mut count = 0;
    for bface in bfaces {
        count += 1;
        writeln!(
            writer,
            "{}\t3\t2\t1\t0\t{}\t{}\t{}\t{}",
            count,
            bface[0] + 1,
            bface[1] + 1,
            bface[2] + 1,
            bface[3] + 1
        )
        .unwrap();
    }
    for (_, nodes) in cubes {
        count += 1;
        writeln!(
            writer,
            "{}\t5\t2\t1\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            count,
            nodes[0] + 1,
            nodes[1] + 1,
            nodes[2] + 1,
            nodes[3] + 1,
            nodes[4] + 1,
            nodes[5] + 1,
            nodes[6] + 1,
            nodes[7] + 1
        )
        .unwrap();
    }
    writeln!(writer, "$EndElements").unwrap();

    writer.flush().unwrap();
}

pub fn vox2txt(voxels: &Array3<i64>, file: &str) {
    let mut f = File::create(file).unwrap();
    for (ijk, v) in voxels.indexed_iter() {
        if *v == 1 {
            writeln!(f, "{}\t{}\t{}", ijk.0, ijk.1, ijk.2).unwrap();
        }
    }
}
