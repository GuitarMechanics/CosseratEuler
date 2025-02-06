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
        materials = ['SUS304', 'NiTi', 'SiliconeRubber']
        materialIndex = int(input('Select material: SUS304[0], NiTi[1], SiliconeRubber[2]'))
        path = f'Matlab/code/csvfiles/{materials[materialIndex]}/*.csv'

        file_list = glob.glob(path)
        self.cbf_list = []
        for file in file_list:
            self.cbf_list.append(cbf(file))
            self.cbf_list.sort(key = lambda x: x.LRratio)

        self.dimensions = []
        for beam in self.cbf_list:
            if (beam.radius, beam.length, beam.LRratio) not in self.dimensions:
                self.dimensions.append((beam.radius, beam.length, beam.LRratio))
        self.dimensions.sort(key = lambda x:x[2])

    def do_plot_all(self, *arg):
        for dim in self.dimensions:
            tmplist = []
            for beam in self.cbf_list:
                if (beam.radius, beam.length) == (dim[0], dim[1]):
                    tmplist.append(beam)
            tmplist.sort(key = lambda x:x.forceratio)
            plt.figure()
            for beam in tmplist:
                plt.subplot(131)
                plt.plot(beam.x, beam.z, label = str(beam.forceratio))
                plt.title("Configuration")
                plt.subplot(132)
                segment, curv = beam.returnCurvature(normalized = True)
                plt.plot(segment, curv, label=str(beam.forceratio))
                plt.title("Curvature")
                plt.subplot(133)
                segment, ang = beam.getAngles(normalized = False)
                plt.plot(segment, ang, label = str(beam.forceratio))
            plt.legend()
            # plt.axis('equal')
            plt.title(f'Radius {dim[0]}, Length {dim[1]}, LR_ratio = {dim[2]}')
            # plt.ylim(0,0.15)
            plt.tight_layout()
            plt.show()
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

    def do_quit(self, *arg):
        return True
    
if __name__ == '__main__':
    MainShell().cmdloop()