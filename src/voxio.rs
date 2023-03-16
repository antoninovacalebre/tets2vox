use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use ndarray::prelude::*;

pub fn vox2gmsh(voxels: &Array3<i64>, dx: f64, file: &str) {
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
            let neighbour_offset = neighbour_offset
                .iter()
                .map(|x| *x as usize)
                .collect::<Vec<usize>>();
            let neighbour = (
                voxel.0 + neighbour_offset[0],
                voxel.1 + neighbour_offset[1],
                voxel.2 + neighbour_offset[2],
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
                    coord[0] + offsets[i][0],
                    coord[1] + offsets[i][1],
                    coord[2] + offsets[i][2],
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
