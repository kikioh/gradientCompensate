function removePath()

% Spinach path
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/experiments'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/kernel'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/etc'));
rmpath(genpath('/pfs/work7/workspace/scratch/ws4078-hmj2/Spinach_dev/interfaces'));

% optimal control path
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/OC/experiments');
rmpath('/pfs/work7/workspace/scratch/ws4078-hmj2/OC/shapes_B0');

end