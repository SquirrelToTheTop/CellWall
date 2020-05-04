import sys
import open3d

def main():
    
    if( len(sys.argv) < 2):
        print("Error no file given !")
        sys.exit()
    
    print("file to process : ", sys.argv[1])

    FOR = open3d.get_o3d_FOR()

    cloud = open3d.read_point_cloud(sys.argv[1]) # Read the point cloud
    open3d.draw_geometries([FOR,cloud]) # Visualize the point cloud     

if __name__ == "__main__":
    main()
