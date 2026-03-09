import SimpleITK as sitk
import os
import numpy as np


def find_moving_img(registerpath):
    swarp_path = os.path.join(registerpath, "s_warp")
    movingImgName = ""
    if os.path.exists(swarp_path):
        files = os.listdir(swarp_path)
        for file in files:
            if file.startswith("sub-") and file.endswith(".nii.gz"):
                movingImgName = file.split(".")[0]
                movingImgName = movingImgName+"_resampled.nii.gz"

    return movingImgName

def load_transform(path):
    try:
        transform_moving2fixed = sitk.ReadTransform(path)
        print("Succesfully loaded transform: " + path.split("/")[-1])
    except:
        print("error loading the transform file")

    return transform_moving2fixed

def load_sitkimage(path):
    try:
        sitkImg = sitk.ReadImage(path)  # Moving image
        print("Succesfully loaded sitk img: " + path.split("/")[-1])
    except:
        print("error loading the sitk img file" + path.split("/")[-1])

    return sitkImg

def map_coordinates(fixedImg, movingImg, transform_moving2fixed):
    nx, ny, nz = fixedImg.GetSize()

    moving_coordinates = np.empty((nx * ny * nz, 3))
    fixed_coordinates = np.empty((nx * ny * nz, 3))

    i = 0
    for x in range(nx):
        if x % 10 == 0:
            print("Progress: " + str(x) + "/" + str(nx))
        for y in range(ny):
            for z in range(nz):
                idx = [x, y, z]
                fixedpnt = fixedImg.TransformIndexToPhysicalPoint(idx)
                movingpnt = transform_moving2fixed.TransformPoint(fixedpnt)
                movingidx = movingImg.TransformPhysicalPointToIndex(movingpnt)
                fixed_coordinates[i, :] = idx
                moving_coordinates[i, :] = movingidx
                i = i + 1

    return moving_coordinates, fixed_coordinates

root="./"
fixedPath=os.path.join(root, "WHS_SD_rat_atlas_v4_pack","WHS_SD_rat_T2star_v1.01.nii.gz")
fixedImg = load_sitkimage(fixedPath)

fixedPoints_name="fixed_img-indeces.npy"
movingPoints_name="moving_img_resampled25um-indeces.npy"

animal=""
path = os.path.join(root, animal, "mri")
sessions = os.listdir(path)

for session in sessions:
    if session.startswith("ses-"):
        print("Animal: "+animal+" session: "+session)
        sessionpath = os.path.join(path, session)
        registerpath = os.path.join(sessionpath, "registration")
        anatpath = os.path.join(sessionpath, "anat")

        if os.path.exists(registerpath):
            print(registerpath)
            transformPath = os.path.join(registerpath, "s_register",
                                         "output_Composite.h5")
            movingImgFilename = find_moving_img(registerpath)
            if movingImgFilename:
                movingPath=os.path.join(anatpath, movingImgFilename)
                fixedPointsSave = os.path.join(registerpath, fixedPoints_name)
                movingPointsSave = os.path.join(registerpath, movingPoints_name)

                if os.path.exists(movingPointsSave):
                    print("file exists, skipping")

                else:
                    print("Mapped-points file does not exists, executing the mapping")
                    movingImg = load_sitkimage(movingPath)
                    transform_moving2fixed = load_transform(transformPath)

                    moving_coordinates, fixed_coordinates = map_coordinates(fixedImg, movingImg, transform_moving2fixed)
                    np.save(fixedPointsSave, fixed_coordinates)
                    np.save(movingPointsSave, moving_coordinates)
