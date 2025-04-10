function removePath()

% Spinach path
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/experiments'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/kernel'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/etc'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/interfaces'));

% spin locking path
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid');
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/experiments');
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/shapes_Bg');
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/OC/shapes_array');
rmpath('/pfs/work7/workspace//scratch/ws4078-hmj2/OC/shapes_array/backup');

end