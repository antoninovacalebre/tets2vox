pub mod progress_bar;
pub mod voxio;

use clap::{arg, command, ArgAction};
use ndarray::prelude::*;

use progress_bar::{pb_init, pb_update};
use voxio::vox2gmsh;

fn main() {
    let start_time = std::time::Instant::now();

    let matches = command!()
        .help_template(
            "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
            ",
        )
        .arg(
            arg!(-i --input ... "Input file")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            arg!(-o --output ... "Output file")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            arg!(-r --res ... "Resolution")
                .required(true)
                .action(ArgAction::Set),
        )
        .get_matches();

    let res: usize = *matches.get_one::<usize>("res").unwrap();
    let input_file = matches.get_one::<String>("input").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();

    let (nodes, tets) = read_gmsh(input_file);

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

    let (vox, h) = tets2vox(&tets_nodes, res);
    let n_vox = vox.iter().filter(|&&x| x == 1).count();

    println!("");
    println!("Voxel size: {}", h);
    println!("Voxel grid shape: {:?}", vox.shape());
    println!("Filled voxels: {}", n_vox);

    // vox2txt(&vox, "tests/voxels.txt");

    vox2gmsh(&vox, h, output_file);

    let end_time = std::time::Instant::now();
    println!("");
    println!(
        "Elapsed time: {} seconds",
        end_time.duration_since(start_time).as_secs_f64()
    );
}

fn tets2vox(tets: &Array3<f64>, res: usize) -> (Array3<i32>, f64) {
    let mut vox: Array3<i32> = Array3::<i32>::zeros((res, res, res));

    // Get the dimensions of the tetrahedra
    let dim: Array1<f64> = tets
        .fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))
        .fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))
        - tets
            .fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y))
            .fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));

    println!("Dimensions: {:?}", dim);

    let dx = dim.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y)) / (res as f64);
    let h = dx.first().unwrap();

    // h is the voxel size (scalar, since the voxels are cubic)

    let mut old_nbars = pb_init("Converting tetrahedra to voxels | 0");
    let min_point = tets
        .fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y))
        .fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y));
    let min_point = min_point;

    for (c, tet) in tets.outer_iter().enumerate() {
        let tet = tet.to_owned() - min_point.to_owned();

        let tet: Array2<f64> = tet.into_shape((4, 3)).unwrap();
        let min_x =
            (tet.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y))[0] / h).round() as usize;
        let max_x =
            (tet.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))[0] / h).round() as usize;
        let min_y =
            (tet.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y))[1] / h).round() as usize;
        let max_y =
            (tet.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))[1] / h).round() as usize;
        let min_z =
            (tet.fold_axis(Axis(0), f64::INFINITY, |&x, &y| x.min(y))[2] / h).round() as usize;
        let max_z =
            (tet.fold_axis(Axis(0), f64::NEG_INFINITY, |&x, &y| x.max(y))[2] / h).round() as usize;

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
