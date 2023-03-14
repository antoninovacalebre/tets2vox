pub mod progress_bar;

use progress_bar::{pb_init, pb_update};
use ndarray::prelude::*;

fn main() {

    let start_time = std::time::Instant::now();

    let (nodes, tets) = read_gmsh("mesh.msh");

    let mut tets_nodes : Array3<f64>;
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

    println!("");
    println!("Voxel size: {}", h);
    println!("Voxel grid shape: {:?}", vox.shape());

    let end_time = std::time::Instant::now();
    println!("");
    println!("Elapsed time: {} seconds", end_time.duration_since(start_time).as_secs_f64());
}

// TODO: implement function vox2gmsh
// fn vox2gmsh(voxels: &Array3<i32>, dx: f64, file: &str){

//     let mappa_r = vec![
//         vec![0, 0, 1, -1, -1, 1, 0, -1, -1, 5, 4],
//         vec![0, 0, -1, 3, 2, -1, -1, 7, 6, -1, -1],
//         vec![0, -1, 0, -1, -1, -1, -1, 0, 1, 2, 3],
//         vec![0, 1, 0, 4, 5, 6, 7, -1, -1, -1, -1],
//         vec![1, 0, 0, -1, 0, 3, -1, -1, 4, 7, -1],
//         vec![-1, 0, 0, 1, -1, -1, 2, 5, -1, -1, 6],
//         vec![-1, 1, -1, 6, -1, -1, -1, -1, -1, -1, -1],
//         vec![1, 1, -1, -1, 7, -1, -1, -1, -1, -1, -1],
//         vec![1, 1, 1, -1, -1, 4, -1, -1, -1, -1, -1],
//         vec![-1, 1, 1, -1, -1, -1, 5, -1, -1, -1, -1],
//         vec![-1, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1],
//         vec![1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1],
//         vec![1, -1, 1, -1, -1, -1, -1, -1, -1, 0, -1],
//         vec![-1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1],
//         vec![0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![-1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![-1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![-1, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![1, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//         vec![0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
//     ];

//     let mappa_o = vec![
//         vec![0, 0, 1],
//         vec![0, 0, -1],
//         vec![0, -1, 0],
//         vec![0, 1, 0],
//         vec![1, 0, 0],
//         vec![-1, 0, 0],
//         vec![-1, 1, -1],
//         vec![1, 1, -1],
//         vec![1, 1, 1],
//         vec![-1, 1, 1],
//         vec![-1, -1, -1],
//         vec![1, -1, -1],
//         vec![1, -1, 1],
//         vec![-1, -1, 1],
//         vec![0, 1, -1],
//         vec![1, 1, 0],
//         vec![0, 1, 1],
//         vec![-1, 1, 0],
//         vec![0, -1, -1],
//         vec![1, -1, 0],
//         vec![0, -1, 1],
//         vec![-1, -1, 0],
//         vec![-1, 0, -1],
//         vec![1, 0, -1],
//         vec![1, 0, 1],
//         vec![-1, 0, 1]
//     ];

//     let offsets = vec![
//         vec![-0.5, 0.5, -0.5],
//         vec![0.5, 0.5, -0.5],
//         vec![0.5, 0.5, 0.5],
//         vec![-0.5, 0.5, 0.5],
//         vec![-0.5, -0.5, -0.5],
//         vec![0.5, -0.5, -0.5],
//         vec![0.5, -0.5, 0.5],
//         vec![-0.5, -0.5, 0.5],
//     ];

//     // Python code to convert to Rust. mappa['o'] is mappa_o, mappa['r'] is mappa_r
//     // all_faces = []
//     // nodes = []
//     // cubes = {}

//     // for (i, j, k), vox in np.ndenumerate(voxels):
//     //     if vox == 0:
//     //         continue
//     //     ijk = np.array([i, j, k])
//     //     # move mesh to origin as minimum and find cube
//     //     coord = (ijk + 0.5) * dx

//     //     cube = np.array([-1, -1, -1, -1, -1, -1, -1, -1])
//     //     for m in mappa:
//     //         if (i+m['o'][0], j+m['o'][1], k+m['o'][2]) in cubes:
//     //             for p, id in enumerate(m['r']):
//     //                 if id >= 0:
//     //                     cube[p] = cubes[(i+m['o'][0], j+m['o'][1], k+m['o'][2])][id]

//     //     for p, c in enumerate(cube):
//     //         if c < 0:
//     //             offset = offsets[p]
//     //             point = coord + offset
//     //             nodes.append(point)
//     //             cube[p] = len(nodes) - 1
        
//     //     cubes[i,j,k] = cube.tolist()

//     //     all_faces.append(cube[np.array([4, 5, 6, 7])])
//     //     all_faces.append(cube[np.array([1, 0, 3, 2])])
//     //     all_faces.append(cube[np.array([1, 5, 6, 2])])
//     //     all_faces.append(cube[np.array([0, 4, 7, 3])])
//     //     all_faces.append(cube[np.array([7, 6, 2, 3])])
//     //     all_faces.append(cube[np.array([0, 1, 5, 4])])

//     let mut all_faces = Vec::new();
//     let mut nodes: Vec<i = Vec::new();
//     let mut cubes = HashMap::new();

//     println!("{}", nodes)

//     for (i, j) in voxels.indexed_iter() {
//         if *j == 0 {
//             continue;
//         }
//         let ijk: Array1<f64> = array![i.0 as f64, i.1 as f64, i.2 as f64];
//         const ONES : Array1<f64> = array![1., 1., 1.];
//         let coord = (ijk + 0.5 * ONES) * dx;

//         let mut cube = array![-1 as i32, -1, -1, -1, -1, -1, -1, -1];
//         for m in mappa_o.iter().zip(mappa_r.iter()) {
//             if let Some(c) = cubes.get(&(i.0 + m.0[0], i.1 + m.0[1], i.2 + m.0[2])) {
//                 for (p, id) in m.1.iter().enumerate() {
//                     if *id >= 0 {
//                         cube[p] = c[*id as usize];
//                     }
//                 }
//             }
//         }

//         for (p, c) in cube.iter_mut().enumerate() {
//             if *c < 0 {
//                 let offset = offsets[p];
//                 let point = coord + offset;
//                 nodes.push(point);
//                 *c = nodes.len() - 1;
//             }
//         }

//         cubes.insert((i.0, i.1, i.2), cube);

//         all_faces.push([cube[4], cube[5], cube[6], cube[7]]);
//         all_faces.push([cube[1], cube[0], cube[3], cube[2]]);
//         all_faces.push([cube[1], cube[5], cube[6], cube[2]]);
//         all_faces.push([cube[0], cube[4], cube[7], cube[3]]);
//         all_faces.push([cube[7], cube[6], cube[2], cube[3]]);
//         all_faces.push([cube[0], cube[1], cube[5], cube[4]]);
//     }

// }

fn tets2vox(tets: &Array3<f64>, res: f64) -> (Array3<i32>, f64) {
    let mut vox: Array3<i32> = Array3::<i32>::zeros((res as usize, res as usize, res as usize));

    let dim = tets.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y)) - tets.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));

    // h is the voxel size (scalar, since the voxels are cubic)
    let dx = dim.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y)) / res;
    let h = dx.first().unwrap();

    let mut old_nbars = pb_init("Converting tetrahedra to voxels | 0");

    for (c, tet) in tets.outer_iter().enumerate() {
        let min_point = tet.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));
        let tet = tet.to_owned() - min_point;
        let tet: Array2<f64> = tet.into_shape((4, 3)).unwrap();
        let min_x = ((tet[[0, 0]].min(tet[[1, 0]]).min(tet[[2, 0]]).min(tet[[3, 0]])) / h).round() as usize;
        let max_x = ((tet[[0, 0]].max(tet[[1, 0]]).max(tet[[2, 0]]).max(tet[[3, 0]])) / h).round() as usize;
        let min_y = ((tet[[0, 1]].min(tet[[1, 1]]).min(tet[[2, 1]]).min(tet[[3, 1]])) / h).round() as usize;
        let max_y = ((tet[[0, 1]].max(tet[[1, 1]]).max(tet[[2, 1]]).max(tet[[3, 1]])) / h).round() as usize;
        let min_z = ((tet[[0, 2]].min(tet[[1, 2]]).min(tet[[2, 2]]).min(tet[[3, 2]])) / h).round() as usize;
        let max_z = ((tet[[0, 2]].max(tet[[1, 2]]).max(tet[[2, 2]]).max(tet[[3, 2]])) / h).round() as usize;

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

        old_nbars = pb_update((c + 1).try_into().unwrap(), tets.shape()[0].try_into().unwrap(), old_nbars, 40);
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

    let out_tets = Array2::from_shape_vec((tets.len(), 4), tets.into_iter().flatten().collect()).unwrap();

    (nodes, out_tets)
}

fn is_point_in_tet(p: &Array1<f64>, tet: &Array2<f64>) -> bool {
    let v = signed_volume(&tet.row(0).to_owned(), &tet.row(1).to_owned(), &tet.row(2).to_owned(), &tet.row(3).to_owned());

    let alpha = signed_volume(&p, &tet.row(1).to_owned(), &tet.row(2).to_owned(), &tet.row(3).to_owned()) / v;
    let beta = signed_volume(&tet.row(0).to_owned(), &p, &tet.row(2).to_owned(), &tet.row(3).to_owned()) / v;
    let gamma = signed_volume(&tet.row(0).to_owned(), &tet.row(1).to_owned(), p, &tet.row(3).to_owned()) / v;
    let delta = signed_volume(&tet.row(0).to_owned(), &tet.row(1).to_owned(), &tet.row(2).to_owned(), &p) / v;

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
