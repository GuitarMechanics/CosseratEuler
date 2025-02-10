from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import glob
import os
import cmd

class MainShell(cmd.Cmd):
    def __init__(self):
        super().__init__()
        os.system('cls' if os.name == 'nt' else 'clear')
        self.intro = 'Cosserat Analyzer'
        self.prompt = '(CosseratAnalyzer)'
        materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']
        materialIndex = int(input('Select material: SUS304[0], NiTi[1], SiliconeRubber[2]'))
        path = f'Matlab/code/csvfiles/{materials[materialIndex]}/*.csv'
        self.materialName = 'NiTi' if materialIndex == 1 else materials[materialIndex]

        file_list = glob.glob(path)
        self.cbf_list = []
        for file in file_list:
            self.cbf_list.append(cbf(file))
            self.cbf_list.sort(key = lambda x: x.LRratio)

        self.dimensions = []
        for beam in self.cbf_list:
            if (beam.radius, beam.length, beam.LRratio) not in self.dimensions:
                self.dimensions.append((beam.radius, beam.length, beam.LRratio))
                print(f'Loaded beam data: L {beam.length} r {beam.radius} LRratio {beam.LRratio}')
        self.dimensions.sort(key = lambda x:x[2])

    def do_plot_all(self, *arg):
        for dim in self.dimensions:
            tmplist = []
            forceratio = 0
            for beam in self.cbf_list:
                if (beam.radius, beam.length) == (dim[0], dim[1]):
                    tmplist.append(beam)
            tmplist.sort(key = lambda x:x.forceratio)
            plt.figure(figsize=(15,6))
            for beam in tmplist:
                forceratio = beam.forceratio
                plt.suptitle(f'{self.materialName}, L = {beam.length} r = {beam.radius} LR_ratio = {beam.LRratio}')
                plt.subplot(131)
                plt.plot(beam.x, beam.z, label = str(beam.forceratio))
                plt.title("Configuration")
                plt.xlabel('x-coordinate')
                plt.ylabel('y-coordinate')
                plt.legend(title='Force ratio')
                plt.subplot(132)
                segment, curv = beam.returnCurvature(normalized = False)
                plt.plot(segment, curv, label=str(beam.forceratio))
                plt.title("Curvature")
                plt.xlabel('Curve length along the track')
                plt.ylabel('Curvature')
                plt.legend(title='Force ratio')
                plt.subplot(133)
                segment, ang = beam.getAngles(normalized = False)
                plt.plot(segment, ang, label = str(beam.forceratio))
                plt.title("Angle of tangent")
                plt.xlabel('Curve length along the track')
                plt.ylabel('Angle of tangent')
                plt.legend(title='Force ratio')
            # plt.axis('equal')
            # plt.ylim(0,0.15)
            plt.tight_layout()
            # plt.show()
            plt.savefig(f'Python/figs/totalfigs/{self.materialName}_r{dim[0]}_l{dim[1]}.png')
            # dummy = input('Draw next?')
    
    def do_plot_curvature_all(self, *arg):
        for dim in self.dimensions:
            tmplist = []
            for beam in self.cbf_list:
                if (beam.radius, beam.length) == (dim[0], dim[1]):
                    tmplist.append(beam)
                tmplist.sort(key = lambda x:x.forceratio)
                plt.figure()
                for beam in tmplist:
                    segment, curv = beam.returnCurvature(normalized = True)
                    plt.plot(segment, curv, label=str(beam.forceratio))
                    plt.title("Curvature")
                plt.legend()
                plt.tight_layout()
                plt.show()
    
    def do_plot_curvature_selected(self, *arg):
        tmplist = []
        radius = float(input('Radius = '))
        length = float(input('length = '))
        normalized = True if input('Normalize?[y/n]') == 'y' else False
        for beam in self.cbf_list:
            if(beam.radius, beam.length) == (radius, length):
                tmplist.append(beam)
        if len(tmplist) == 0:
            print("No beams in such dimensions\n")
            return
        tmplist.sort(key = lambda x : x.forceratio)
        plt.figure()
        for beam in tmplist:
            segment, curv = beam.returnCurvature(normalized)
            plt.plot(segment, curv, label = str(beam.forceratio))
        plt.title(f'Curvature: r = {radius}, L = {length}')
        plt.xlabel('Position on the curve')
        plt.ylabel('Curvature')
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    def do_plot_configuration_selected(self, *arg):
        tmplist = []
        radius = float(input('Radius = '))
        length = float(input('length = '))
        normalized = True if input('Normalize?[y/n]') == 'y' else False
        for beam in self.cbf_list:
            if(beam.radius, beam.length) == (radius, length):
                tmplist.append(beam)
        if len(tmplist) == 0:
            print("No beams in such dimensions\n")
            return
        tmplist.sort(key = lambda x : x.forceratio)
        plt.figure()
        for beam in tmplist:
            x = beam.x / beam.length if normalized else beam.x
            y = beam.z / beam.length if normalized else beam.z
            plt.plot(x, y, label = str(beam.forceratio))
        plt.title(f'Configuration: r = {radius}, L = {length}, Normalized: {normalized}')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def do_plot_angle_selected(self, *arg):
        tmplist = []
        radius = float(input('Radius = '))
        length = float(input('length = '))
        normalized = True if input('Normalize?[y/n]') == 'y' else False
        for beam in self.cbf_list:
            if(beam.radius, beam.length) == (radius, length):
                tmplist.append(beam)
        if len(tmplist) == 0:
            print("No beams in such dimensions\n")
            return
        tmplist.sort(key = lambda x : x.forceratio)
        plt.figure()
        for beam in tmplist:
            segment, curv = beam.getAngles(normalized)
            plt.plot(segment, curv, label = str(beam.forceratio))
        plt.title(f'Angle: r = {radius}, L = {length}, Normalize: {normalized}')
        plt.xlabel('Position on the curve')
        plt.ylabel('Angle(deg)')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def do_plot_normalized_curvature_all(self, *arg):
        lr_ratio_upper = float(input('Max. l_r ratio? \n'))
        lr_ratio_lower = float(input('Min. l_r ratio? \n'))
        
        plt.figure()
        for beam in self.cbf_list:
            if beam.LRratio <= lr_ratio_upper and beam.LRratio >= lr_ratio_lower:
                x, y = beam.getCurvature(normalize = True)
                print(f'plotted curvature: L {beam.length} r {beam.radius} Fr {beam.forceratio}')
                plt.plot(x,y)
        plt.title(f'Curvature plot: LR_ratio {lr_ratio_lower} to {lr_ratio_upper}')
        plt.show()


    def do_quit(self, *arg):
        return True
    
if __name__ == '__main__':
    MainShell().cmdloop()