pub mod progress_bar;

use ndarray::prelude::*;
use progress_bar::{pb_init, pb_update};
use std::collections::HashMap;
use std::fs::File;

fn main() {
    let start_time = std::time::Instant::now();

    let (nodes, tets) = read_gmsh("tests/mesh.msh");

    let mut tets_nodes: Array3<f64>;
    tets_nodes = Array3::<f64>::zeros((tets.shape()[0], 4, 3));
    for (i, tet) in tets.outer_iter().enumerate() {
        for (j, node) in tet.iter().enumerate() {
            tets_nodes[[i, j, 0]] = nodes[[*node as usize, 0]];
            tets_nodes[[i, j, 1]] = nodes[[*node as usize, 1]];
            tets_nodes[[i, j, 2]] = nodes[[*node as usize, 2]];
        }
    }

    println!("");
    println!("Number of tetrahedra: {}", tets_nodes.shape()[0]);
    println!("Number of nodes: {}", nodes.shape()[0]);

    let (vox, h) = tets2vox(&tets_nodes, 100.0);
    let n_vox = vox.iter().filter(|&&x| x == 1).count();

    println!("");
    println!("Voxel size: {}", h);
    println!("Voxel grid shape: {:?}", vox.shape());
    println!("Filled voxels: {}", n_vox);

    let end_time = std::time::Instant::now();
    println!("");
    println!(
        "Elapsed time: {} seconds",
        end_time.duration_since(start_time).as_secs_f64()
    );
}

fn vox2gmsh(voxels: &Array3<i32>, dx: f64, file: &str) {
    // Constants

    // Map of faces to stitch together as a list
    // of tuples (vector_offsets, faces_to_stitch)
    let stitch_map: Vec<(Vec<i32>, Vec<i32>)> = vec![
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

    let ncell = voxels.iter().filter(|&&x| x != 0).count();

    let mut all_faces: Vec<Vec<i32>> = Vec::with_capacity(ncell * 6);
    let mut nodes = vec![];
    let mut cubes: HashMap<(usize, usize, usize), Vec<i32>> = HashMap::new();

    for (ijk, v) in voxels.indexed_iter() {
        if *v == 0 {
            continue;
        }

        // move mesh to origin as minimum and find cube
        let mut coord = vec![ijk.0 as f64, ijk.1 as f64, ijk.2 as f64];
        coord = coord.iter().map(|x| (x + 0.5) * dx).collect();

        // Initialize cube as 1D vector of 8 nodes (-1)
        let mut cube = Array::from_elem((8,), -1);

        // Iterate over stich_map
        for (_, (stitch, faces)) in stitch_map.iter().enumerate() {
            // check if tuple is a key in cubes
            let stitch = stitch.iter().map(|x| *x as usize).collect::<Vec<usize>>();
            let tup = (ijk.0 + stitch[0], ijk.1 + stitch[1], ijk.2 + stitch[2]);
            if cubes.contains_key(&ijk) {
                // iterate over faces
                for (j, face) in faces.iter().enumerate() {
                    if *face >= 0 {
                        cube[j] = cubes[&tup][*face as usize];
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
                cube[i] = nodes.len() as i32 - 1;
            }
        }

        cubes.insert(ijk, cube.to_vec());

        all_faces.push(vec![cube[4], cube[5], cube[6], cube[7]]);
        all_faces.push(vec![cube[1], cube[0], cube[3], cube[2]]);
        all_faces.push(vec![cube[1], cube[5], cube[6], cube[2]]);
        all_faces.push(vec![cube[0], cube[4], cube[7], cube[3]]);
        all_faces.push(vec![cube[7], cube[6], cube[2], cube[3]]);
        all_faces.push(vec![cube[0], cube[1], cube[5], cube[4]]);
    }

    // Filter all faces that have no duplicates
    let mut bfaces = all_faces
        .into_iter() // convert to iterator
        .map(|mut face| {
            face.sort();
            face
        }) // sort the elements of each face
        .collect::<Vec<Vec<i32>>>()
        .into_iter()
        .collect::<std::collections::HashSet<Vec<i32>>>() // convert to hash set
        .into_iter()
        .collect::<Vec<Vec<i32>>>(); // convert back to vector

    // # export mesh
    // with open(file, 'w') as f:
    //     f.write('$MeshFormat\n')
    //     f.write('2.2 0 8\n')
    //     f.write('$EndMeshFormat\n')
    //     f.write('$Nodes\n')
    //     f.write(f'{len(nodes)}\n')
    //     for i, mnode in enumerate(nodes):
    //         f.write(f'{i+1}\t{mnode[0]}\t{mnode[1]}\t{mnode[2]}\n')
    //     f.write('$EndNodes\n')
    //     f.write('$Elements\n')
    //     f.write(f'{len(cubes)+len(bfaces)}\n')
    //     count = 0
    //     for bface in bfaces:
    //         count += 1
    //         f.write(f'{count}\t3\t2\t1\t0\t{bface[0]+1}\t{bface[1]+1}\t{bface[2]+1}\t{bface[3]+1}\n')
    //     for xyz, nodes in cubes.items():
    //         count += 1
    //         f.write(f'{count}\t5\t2\t1\t0\t{nodes[0]+1}\t{nodes[1]+1}\t{nodes[2]+1}\t{nodes[3]+1}\t{nodes[4]+1}\t{nodes[5]+1}\t{nodes[6]+1}\t{nodes[7]+1}\n')
    //     f.write('$EndElements\n')
    //     f.write('\n')

    // open file
    let mut f = File::create(file).unwrap();
    



}

fn tets2vox(tets: &Array3<f64>, res: f64) -> (Array3<i32>, f64) {
    let mut vox: Array3<i32> = Array3::<i32>::zeros((res as usize, res as usize, res as usize));

    let dim = tets.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))
        - tets.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));

    // h is the voxel size (scalar, since the voxels are cubic)
    let dx = dim.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y)) / res;
    let h = dx.first().unwrap();

    let mut old_nbars = pb_init("Converting tetrahedra to voxels | 0");

    for (c, tet) in tets.outer_iter().enumerate() {
        let min_point = tet.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));
        let tet = tet.to_owned() - min_point;
        let tet: Array2<f64> = tet.into_shape((4, 3)).unwrap();
        let min_x = ((tet[[0, 0]]
            .min(tet[[1, 0]])
            .min(tet[[2, 0]])
            .min(tet[[3, 0]]))
            / h)
            .round() as usize;
        let max_x = ((tet[[0, 0]]
            .max(tet[[1, 0]])
            .max(tet[[2, 0]])
            .max(tet[[3, 0]]))
            / h)
            .round() as usize;
        let min_y = ((tet[[0, 1]]
            .min(tet[[1, 1]])
            .min(tet[[2, 1]])
            .min(tet[[3, 1]]))
            / h)
            .round() as usize;
        let max_y = ((tet[[0, 1]]
            .max(tet[[1, 1]])
            .max(tet[[2, 1]])
            .max(tet[[3, 1]]))
            / h)
            .round() as usize;
        let min_z = ((tet[[0, 2]]
            .min(tet[[1, 2]])
            .min(tet[[2, 2]])
            .min(tet[[3, 2]]))
            / h)
            .round() as usize;
        let max_z = ((tet[[0, 2]]
            .max(tet[[1, 2]])
            .max(tet[[2, 2]])
            .max(tet[[3, 2]]))
            / h)
            .round() as usize;

        for i in min_x..max_x {
            for j in min_y..max_y {
                for k in min_z..max_z {
                    if vox[[i, j, k]] == 1 {
                        continue;
                    }

                    let x = (i as f64 + 0.5) * h;
                    let y = (j as f64 + 0.5) * h;
                    let z = (k as f64 + 0.5) * h;

                    let point = array![x, y, z];

                    if is_point_in_tet(&point, &tet) {
                        vox[[i, j, k]] = 1;
                    }
                }
            }
        }

        old_nbars = pb_update(
            (c + 1).try_into().unwrap(),
            tets.shape()[0].try_into().unwrap(),
            old_nbars,
            40,
        );
    }

    (vox, *h)
}

fn read_gmsh(file: &str) -> (Array2<f64>, Array2<i32>) {
    let lines = std::fs::read_to_string(file).unwrap();
    let lines: Vec<&str> = lines.split("\r").collect();

    if lines[1].split(" ").collect::<Vec<&str>>()[0].trim() != "2.2" {
        panic!("Only Gmsh 2.2 is supported");
    }

    let nnode = lines[4].trim().parse::<usize>().unwrap();
    let mut nodes = Array2::<f64>::zeros((nnode, 3));

    for i in 0..nnode {
        let line = lines[5 + i].split(" ").collect::<Vec<&str>>();
        nodes[[i, 0]] = line[1].trim().parse::<f64>().unwrap();
        nodes[[i, 1]] = line[2].trim().parse::<f64>().unwrap();
        nodes[[i, 2]] = line[3].trim().parse::<f64>().unwrap();
    }

    let nelem = lines[7 + nnode].trim().parse::<usize>().unwrap();
    let mut tets = Vec::new();

    for i in 0..nelem {
        let line = lines[8 + nnode + i].split(" ").collect::<Vec<&str>>();
        if line[1] == "4" {
            tets.push([
                line[5].trim().parse::<i32>().unwrap() - 1,
                line[6].trim().parse::<i32>().unwrap() - 1,
                line[7].trim().parse::<i32>().unwrap() - 1,
                line[8].trim().parse::<i32>().unwrap() - 1,
            ]);
        }
    }

    let out_tets =
        Array2::from_shape_vec((tets.len(), 4), tets.into_iter().flatten().collect()).unwrap();

    (nodes, out_tets)
}

fn is_point_in_tet(p: &Array1<f64>, tet: &Array2<f64>) -> bool {
    let v = signed_volume(
        &tet.row(0).to_owned(),
        &tet.row(1).to_owned(),
        &tet.row(2).to_owned(),
        &tet.row(3).to_owned(),
    );

    let alpha = signed_volume(
        &p,
        &tet.row(1).to_owned(),
        &tet.row(2).to_owned(),
        &tet.row(3).to_owned(),
    ) / v;
    let beta = signed_volume(
        &tet.row(0).to_owned(),
        &p,
        &tet.row(2).to_owned(),
        &tet.row(3).to_owned(),
    ) / v;
    let gamma = signed_volume(
        &tet.row(0).to_owned(),
        &tet.row(1).to_owned(),
        p,
        &tet.row(3).to_owned(),
    ) / v;
    let delta = signed_volume(
        &tet.row(0).to_owned(),
        &tet.row(1).to_owned(),
        &tet.row(2).to_owned(),
        &p,
    ) / v;

    alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0 && delta >= 0.0
}

fn signed_volume(a: &Array1<f64>, b: &Array1<f64>, c: &Array1<f64>, d: &Array1<f64>) -> f64 {
    let ab = b - d;
    let ac = c - d;
    let ad = a - d;
    let cross = cross_product(&ab, &ac);
    let volume = ad.dot(&cross) / 6.0;
    volume
}

fn cross_product(a: &Array1<f64>, b: &Array1<f64>) -> Array1<f64> {
    if a.len() != 3 || b.len() != 3 {
        panic!("Cross product only defined for 3D vectors");
    }

    let mut cross = Array1::<f64>::zeros(3);
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
    cross
}
