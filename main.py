# from PyQt5 import QtGui
from PyQt5.QtWidgets import QMainWindow,QApplication,QMessageBox , QVBoxLayout , QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.pyplot  import Figure
from TraverseCurve import pressure_traverse
from ipr import vogel_fe_undersaturated_fes,vogel_fe_saturated_fes,process_vlp
# import matplotlib.pyplot as plt
from PyQt5.uic import loadUi
from sys import argv
# import numpy as np
from numpy import arange, linspace
# ------------------------- the canvas for the graphs ---------------------

class MatplotlibWidget(QWidget):
    def __init__(self,parent =None):
            super(MatplotlibWidget,self).__init__(parent)
            self.figure = Figure()
            self.canvas = FigureCanvasQTAgg(self.figure)
            self.axis=self.figure.add_subplot(111)
            self.layoutvertical= QVBoxLayout(self)
            self.layoutvertical.addWidget(self.canvas)
            # initiate the grapfh
            self.figure.gca().invert_yaxis()
            # self.axis.xaxis.tick_top()
            # self.axis.xaxis.set_label_position("top")
            self.axis.set_xlabel("Pressure (psi)")
            self.axis.set_ylabel("Depth (ft)")
            
class VLPWidget(QWidget):
    def __init__(self,parent =None):
            super(VLPWidget,self).__init__(parent)
            self.figure = Figure()
            self.canvas = FigureCanvasQTAgg(self.figure)
            self.axis=self.figure.add_subplot(111)
            self.layoutvertical= QVBoxLayout(self)
            self.layoutvertical.addWidget(self.canvas)
            # initiate the grapfh
            # self.figure.gca().invert_yaxis()
            # self.axis.xaxis.tick_top()
            # self.axis.xaxis.set_label_position("top")
            self.axis.set_xlabel("Q (STB/Day)")
            self.axis.set_ylabel("Pwf (psi)")
            
# ----------- main class -------------------------
class MainApp(QMainWindow):
    def __init__(self):
        super().__init__()
        # initiate the gui 
        loadUi("Ui/traverse_curve.ui",self)
        # intiate the ui setting
        self.ui_setting()
        self.buttons()
        # self.setup_icons()
        self.init_widget()
        
    def ui_setting(self):
        self.setWindowTitle("Pressure Traverse Crve")
        # init the drawing area 
        
    def init_widget(self):
        # for the traverse curve
        self.tc = MatplotlibWidget()
        self.layoutvertical= QVBoxLayout(self.widget_2)
        self.layoutvertical.addWidget(self.tc)
        # for the vlp
        self.vlp = VLPWidget()
        self.layoutvertical= QVBoxLayout(self.widget)
        self.layoutvertical.addWidget(self.vlp)
        
        # IPR_setting
        self.ipr = VLPWidget()
        self.layoutvertical= QVBoxLayout(self.widget_3)
        self.layoutvertical.addWidget(self.ipr)
        
        # nodal analysis tab widget
        self.nodal_a = VLPWidget()
        self.layoutvertical= QVBoxLayout(self.widget_4)
        self.layoutvertical.addWidget(self.nodal_a)
        
    def buttons(self):
        self.pushButton_3.clicked.connect(self.plot_curve)
        self.pushButton_2.clicked.connect(self.show_multiple_curves)
        self.pushButton.clicked.connect(self.show_vlp)
        self.pushButton_4.clicked.connect(self.show_ipr)
        self.pushButton_5.clicked.connect(self.nodal_point)
        
    def inputs(self):
        try:
            self.liquid_rate = float(self.l7.text())
            self.wc = float(self.l12.text())
            self.GLR = float(self.l8.text())
            self.gas_grav = float(self.l10.text())
            self.oil_grav = float(self.l9.text())
            self.wtr_grav = float(self.l11.text())
            self.diameter = float(self.l5.text())
            self.angle = float(self.l13.text())
            self.margin = 50
            self.thp = float(self.l2.text())
            self.tht = float(self.l3.text())
            self.twf = float(self.l4.text())
            self.total_depth = float(self.l1.text())
            self.sample_size = 100
        except:
            QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
    
    def inputs_ipr(self):
        try:
            pr = float(self.l1_2.text())
            pb = float(self.l2_2.text())
            fe = float(self.l3_2.text())
            pwf = float(self.l4_2.text())
            q = float(self.l5_2.text())
            
            return [pr,pb,fe,pwf,q]
        except:
            QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
        
    def plot_curve(self,rate=False,G=False,Wc=False):
        try:
            self.tc.figure.legend("",frameon=False)
            self.inputs()
            self.tc.axis.clear()
            p, depths = pressure_traverse(liquid_rate=self.liquid_rate, depth=self.total_depth, tht=self.tht, twf=self.twf,
                                            glr=self.GLR, wc=self.wc, gas_grav=self.gas_grav, oil_grav=self.oil_grav,
                                            wtr_grav=self.wtr_grav, diameter=self.diameter, angle=self.angle, margin=self.margin,
                                            thp=self.thp, sample_size=self.sample_size)
            if rate:
                self.tc.axis.plot(p, depths, label=f'Liquid Rate = {self.liquid_rate} STB\day \n')
            elif G:
                self.tc.axis.plot(p, depths, label=f'GLR =  {self.GLR} scf\stb \n')
            elif Wc:
                self.tc.axis.plot(p, depths, label=f'WC =  {self.GLR} ')

            else:
                self.tc.axis.plot(p, depths, label=f'D = {self.diameter} inch \n'
                                        f'Liquid Rate = {self.liquid_rate} STB\day \n'
                                        f'GLR =  {self.GLR} scf\stb \n'
                                        f'wc = {self.wc}\nGas gravity = {self.gas_grav} \n'
                                        f'Oil gravity=  {self.oil_grav} API\n'
                                        f'Water Gravity = {self.wtr_grav}\n\n')
            # Show the minor grid lines with very faint and almost transparent grey lines
            self.tc.axis.minorticks_on()
            self.tc.axis.grid(which='major', color='#666666', linestyle='-')
            self.tc.axis.grid(which='minor', color='#999999', linestyle='-',alpha=0.3)
            self.tc.axis.set_ylim(0, depths[-1])
            self.tc.figure.gca().invert_yaxis()
            self.tc.axis.set_xlabel("Pressure (psi)")
            self.tc.axis.set_ylabel("Depth (ft)")
            self.tc.axis.legend(loc="upper right",framealpha=.5) 
            self.tc.canvas.draw()
        except:
            QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
        
    def show_multiple_curves(self):
        try:
            self.tc.figure.legend("",frameon=False)
            self.tc.axis.clear()
            # chekc which variable it is
            variable = self.comboBox.currentText()
            # get the min and max and points 
            min_value = float(self.l14.text())
            max_value = float(self.l15.text())
            step = float(self.l16.text())
            values = list(arange(min_value, max_value ,step))+[max_value]
            
            if variable  == "Multitple Flow rates":
                for q in values :
                    p, depths = pressure_traverse(liquid_rate=q, depth=self.total_depth, tht=self.tht, twf=self.twf,
                                            glr=self.GLR, wc=self.wc, gas_grav=self.gas_grav, oil_grav=self.oil_grav,
                                            wtr_grav=self.wtr_grav, diameter=self.diameter, angle=self.angle, margin=self.margin,
                                            thp=self.thp, sample_size=self.sample_size)
                    self.tc.axis.plot(p, depths, label=f'Liquid Rate = {round(q,2)} STB\day ')
            
            elif variable == "Multiple GLRs":
                for glr in values :
                    p, depths = pressure_traverse(liquid_rate=self.liquid_rate, depth=self.total_depth, tht=self.tht, twf=self.twf,
                                            glr=glr, wc=self.wc, gas_grav=self.gas_grav, oil_grav=self.oil_grav,
                                            wtr_grav=self.wtr_grav, diameter=self.diameter, angle=self.angle, margin=self.margin,
                                            thp=self.thp, sample_size=self.sample_size)
                    self.tc.axis.plot(p, depths, label=f'GLR =  {round(glr,2)} scf\stb ')
                
            elif variable == "Multiple WC":
                for wc in values :
                    p, depths = pressure_traverse(liquid_rate=self.liquid_rate, depth=self.total_depth, tht=self.tht, twf=self.twf,
                                            glr=self.GLR, wc=wc, gas_grav=self.gas_grav, oil_grav=self.oil_grav,
                                            wtr_grav=self.wtr_grav, diameter=self.diameter, angle=self.angle, margin=self.margin,
                                            thp=self.thp, sample_size=self.sample_size)
                    self.tc.axis.plot(p, depths, label=f'WC =  {round(wc,2)} %')
                
            # Show the minor grid lines with very faint and almost transparent grey lines
            self.tc.axis.minorticks_on()
            self.tc.axis.grid(which='major', color='#666666', linestyle='-')
            self.tc.axis.grid(which='minor', color='#999999', linestyle='-',alpha=0.3)
            self.tc.axis.set_ylim(0, depths[-1])
            self.tc.figure.gca().invert_yaxis()
            self.tc.axis.set_xlabel("Pressure (psi)")
            self.tc.axis.set_ylabel("Depth (ft)")
            self.tc.axis.legend(loc="upper right",framealpha=.5) 
            self.tc.canvas.draw()
        except:
            QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
            
    def show_vlp(self):
        try:
            self.vlp.axis.legend().remove()
            self.inputs()
            self.vlp.axis.clear()
            minq = float(self.l17.text())
            maxq = float(self.l18.text())
            rates = linspace(minq, maxq,100)
            bhps = []
            for q in rates:
                p, _ = pressure_traverse(liquid_rate=q, depth=self.total_depth, tht=self.tht, twf=self.twf,
                                            glr=self.GLR, wc=self.wc, gas_grav=self.gas_grav, oil_grav=self.oil_grav,
                                            wtr_grav=self.wtr_grav, diameter=self.diameter, angle=self.angle, margin=self.margin,
                                            thp=self.thp, sample_size=self.sample_size)
                bhp = p[-1]
                bhps.append(bhp)
            
            # definig the vlp data 
            self.vlp_data = [rates,bhps]
                
            self.vlp.axis.plot(rates, bhps, "-")
            self.vlp.axis.set_xlabel('Flow rate STB/Day')
            self.vlp.axis.set_ylabel('PWf Psi')
            
            # Show the minor grid lines with very faint and almost transparent grey lines
            self.vlp.axis.minorticks_on()
            self.vlp.axis.grid(which='major', color='#666666', linestyle='-')
            self.vlp.axis.grid(which='minor', color='#999999', linestyle='-',alpha=0.5)
            self.vlp.canvas.draw()
        except:
            QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
            
    def show_ipr(self):
        # try:
        # self.ipr.axis.legend().remove()
        pr,pb,fe,pwf,q = self.inputs_ipr()
        self.ipr.axis.clear()
        
        # check if the reservoir is saturated or undersaturated
        if pr < pb:
            Qs,Ps =  vogel_fe_saturated_fes(pr,pwf,q,fe,[fe],20)
        else:
            Qs,Ps = vogel_fe_undersaturated_fes(pr,pb,pwf,q,fe,[fe],20)
        
        Q0 =Qs[0]
        P0 = Ps[0]
        
        # defining the IPR attrebute to the app
        self.ipr_data = [Q0,P0]
            
        self.ipr.axis.plot(Q0, P0, "-")
        self.ipr.axis.set_xlabel('Flow rate STB/Day')
        self.ipr.axis.set_ylabel('PWf Psi')
        
        # Show the minor grid lines with very faint and almost transparent grey lines
        self.ipr.axis.minorticks_on()
        self.ipr.axis.grid(which='major', color='#666666', linestyle='-')
        self.ipr.axis.grid(which='minor', color='#999999', linestyle='-',alpha=0.5)
        self.ipr.canvas.draw()
        # except:
        #     QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
    def nodal_point(self):
        # try:
        self.nodal_a.axis.legend().remove()
        self.nodal_a.axis.clear()
        # checking the ipr and vlp data are defined
        if not (self.ipr_data and self.vlp_data ):
            raise exception("You must first get the data from the IPR and VLP.")
        
        Q_i,P_i = self.ipr_data
        vlp_data = process_vlp(self.vlp_data)
        Q_v,P_v = vlp_data
        self.nodal_a.axis.plot(Q_i, P_i, "-",label="Inflow")
        self.nodal_a.axis.plot(Q_v, P_v, "-",label="outflow")
        self.nodal_a.axis.legend(loc="upper right",framealpha=.5) 
        self.nodal_a.axis.set_xlabel('Flow rate STB/Day')
        self.nodal_a.axis.set_ylabel('PWf Psi')
        
        # Show the minor grid lines with very faint and almost transparent grey lines
        self.nodal_a.axis.minorticks_on()
        self.nodal_a.axis.grid(which='major', color='#666666', linestyle='-')
        self.nodal_a.axis.grid(which='minor', color='#999999', linestyle='-',alpha=0.5)
        self.nodal_a.canvas.draw()
        # except:
        #     QMessageBox.warning(self,"Warning","There is some error, make sure that all variables are entred right.")
        
    
def main():
    app =QApplication(argv)
    window = MainApp()
    window.show()
    app.exec_()
if __name__=="__main__":
    main()