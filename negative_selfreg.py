from manim import *
from bio_v2 import *
import numpy as np
from collections import defaultdict
from statistics import variance, mean

config.background_color = WHITE

class Negative(Scene):

    def construct(self):
        self.beta_vary = [1,3]
        #network = self.simple_regulation()
        #self.play(network.animate.scale(0.5),network.animate.move_to(2*LEFT+2*UP))
        self.graph_scene()
        #self.simple_graph()
        self.play(Create(self.neg_axes),Create(self.neg_labels),
                  Create(self.simple_axes),Create(self.simple_labels),FadeIn(self.title_rate),
                  FadeIn(self.title_neg),FadeIn(self.title_simple))
        cells = self.arrange_cells()
        stable_states = defaultdict(list)
        projections = defaultdict(lambda: defaultdict(list))

        simple_g = self.simple_reg_graph()
        neg_g = self.negative_auto_graph()
        Xst = bisection(lambda x_lvl: self.neg_protein.production_rate(x_lvl),
                        lambda x_lvl: self.neg_protein.degradation_rate(x_lvl),
                        0, 10)
        stable_states['neg'].append(Xst)
        projections['neg'][self.neg_protein.production_rates['max']] = self.get_projection(self.neg_protein)
        Xst = bisection(lambda x_lvl: self.simple_protein.production_rate(x_lvl),
                        lambda x_lvl: self.simple_protein.degradation_rate(x_lvl),
                        0, 10)
        stable_states['simple'].append(Xst)
        projections['simple'][self.simple_protein.production_rates['max']] = self.get_projection(self.simple_protein)

        # plot main graph
        for (k, v), (k2, v2) in zip(simple_g.items(), neg_g.items()):
            if k == 'K':
                continue
            elif k == 'stable':
                self.play(Create(v[0]), Create(v2[0]), Create(v[1]), Create(v2[1]))
            else:
                self.play(Create(v[0]), Create(v2[0]))
                self.wait(0)

        # plot cell variations of protein
        neg_vary_g,neg_vary_p = self.get_variations_protein(self.neg_protein,self.beta_vary)
        simple_vary_g,simple_vary_p = self.get_variations_protein(self.simple_protein,self.beta_vary)
        for beta in self.beta_vary:
            for (k, v), (k2, v2) in zip(list(simple_vary_g[beta].items()),list(neg_vary_g[beta].items())):
                if k == 'K':
                    continue
                elif k == 'stable':
                    self.play(Create(v[0]), Create(v2[0]),Create(v[1]), Create(v2[1]))
                else:
                    self.play(Create(v[0]), Create(v2[0]))
                    self.wait(0)

        # for beta in self.beta_vary:
            protein = neg_vary_p[beta]
            Xst = bisection(lambda x_lvl: protein.production_rate(x_lvl),
                            lambda x_lvl: protein.degradation_rate(x_lvl),
                            0, 10)
            stable_states['neg'].append(Xst)
            projections['neg'][beta] = self.get_projection(protein)
            protein = simple_vary_p[beta]
            Xst = bisection(lambda x_lvl: protein.production_rate(x_lvl),
                            lambda x_lvl: protein.degradation_rate(x_lvl),
                            0, 10)
            stable_states['simple'].append(Xst)
            projections['simple'][beta] = self.get_projection(protein)

        # label things
        b1 = self.brace_label(list(x[1] for x in projections['neg'].values()),
                              "stable states in neg. autoreg.",DOWN)
        b2 = self.brace_label(list(x[1] for x in projections['simple'].values()),
                              "stable states in simple reg.",DOWN)
        b3 = self.brace_label([x[0] for beta in self.beta_vary for x in neg_vary_g[beta].values()],
                         "production levels",RIGHT)
        b4 = self.brace_label([x[0] for beta in self.beta_vary for x in simple_vary_g[beta].values()],
                         "production levels", RIGHT)
        self.Appear(VGroup(*[b1,b2,b3,b4]))
        #plot distribution of variations
        #self.plot_normal_graph(stable_states,projections)
        self.wait()

    def graph_scene(self):
        self.number_plane = NumberPlane(
            x_range=[-10, 10, 1],
            y_range=[-10, 10, 1],
            background_line_style={
                "stroke_color": BLACK,
                "stroke_width": 1,
                "stroke_opacity": 0.4
            })
        self.number_plane.coords_to_point(0)
        # rate of change plots
        n = self.get_axes()
        s = self.get_axes()
        self.neg_axes, self.neg_labels = n[0], n[1]
        self.simple_axes,self.simple_labels = s[0],s[1]

        # distribution plots
        normal_n = self.get_axes(x_label="steady/states",
                                 y_label="probability",
                                 xr=[-10,10,1],yr=[0,1,0.1],scale=0.35,color=WHITE,ticks=False)
        normal_s = self.get_axes(x_label="steady/states",
                                 y_label="probability",
                                 xr=[-10,10,1],yr=[0,1,0.1],scale=0.35,color=WHITE,ticks=False)
        self.nnormal_axes,self.nnormal_labels = normal_n[0],normal_n[1]
        self.snormal_axes,self.snormal_labels = normal_s[0],normal_s[1]

        # arrange in the grid
        self.graphs_group = Group(s,normal_s,n,normal_n).arrange_in_grid(2,2,buff=(0.5,0.5))

        #graph titles
        self.title_rate = Tex("Rate of change",
                         color=BLACK).scale(0.65).next_to(self.simple_axes.get_top()+LEFT)
        self.title_dist = Tex("Distribution of stable states",
                         color=BLACK).scale(0.65).move_to(RIGHT*2.5).align_to(self.title_rate,direction=UP)
        self.title_simple =  Tex("Simple",
                            color=BLACK).scale(0.5).rotate(PI/2.3).next_to(self.simple_axes.get_right(),buff=0.5)
        self.title_neg = Tex("Negative \n autoregulation",
                        color=BLACK).rotate(PI/2.3).scale(0.5).next_to(self.neg_axes.get_right(),buff=0.5)


    def get_variations_protein(self,protein,variations=[1,6]):
        # return list of naturally occuring variations of a given protein
        colors = [PURPLE,PINK]
        if type(protein) == ProteinNegativeAuto:
            vary_plots = {}
            vary_proteins = {}
            axes = self.neg_axes
            for i,beta in enumerate(variations):
                p = ProteinNegativeAuto(start=0, leaky_production_rate=0,
                                        max_production_rate=beta,
                                        degradation_rate=protein.alpha)
                p.add_repressor(p, protein.r_kinetics['K'], protein.r_kinetics['n'])
                plot = self.get_plots(p, axes=axes, plot_types=set(['prod', 'max', 'steady']),col=colors[i])
                vary_plots[beta] = plot
                vary_proteins[beta] = p
        elif type(protein) == ProteinSimpleReg:
            vary_plots = {}
            vary_proteins = {}
            axes = self.simple_axes
            for i,beta in enumerate(variations):
                p = ProteinSimpleReg(start=0,max_production_rate=beta,
                                        degradation_rate=protein.alpha)
                plot = self.get_plots(p, axes=axes, plot_types=set(['prod', 'max', 'steady']),col=colors[i])
                vary_plots[beta] = plot
                vary_proteins[beta]= p
        return vary_plots, vary_proteins

    def get_axes(self,x_label="Y",y_label=r"\frac{dY}{dt}",xr=[0,7,1],yr=[0,7,1],color=BLACK,scale=0.5,ticks=True):
        axes = Axes(
            x_range=xr,
            y_range=yr,
            axis_config={"color": color, "include_tip": False,"include_ticks": ticks},
        ).scale(scale)
        axes.coords_to_point(0)
        y_label = axes.get_y_axis_label(MathTex(y_label), edge=LEFT, direction=LEFT, buff=0.4).set_color(
            BLACK).scale(scale)
        x_label = axes.get_x_axis_label(x_label).set_color(BLACK).scale(scale)
        labels = VGroup(x_label, y_label)
        return Group(axes,labels)

    def plot_normal_graph(self,stable_states=None,projections=None):
        """
        TODO: show relationship between variance here and  on protein plots
        1. when varying production, get projection
        2. transform projections into distribution
        3. show variance on both distributions
        4. one bigger than the other!

        :param stable_states: {'neg':[],'simp':[]}
        :param projections: {'neg':[],'simp':[]}
        :return:
        """
        # Creating a series of data of in range of 1-50.
        def gauss(x,mid=0,variance=1, h=1):
            # compute mid and variance from stable states
            from math import pow, exp
            return h * exp(-pow(x - mid, 2) / (2 * variance))

        # get pdf and transformations of projections for negative autoregulation
        mean_n,var_n = mean(stable_states['neg']), variance(stable_states['neg'])
        pdf_n = self.nnormal_axes.get_graph(lambda x_lvl: gauss(x_lvl,mean_n,var_n),
                       color=RED,
                       x_range=[-2, 8])
        proj_trans_n = []
        for x in stable_states['neg']:
            y = gauss(x,mean_n,var_n)
            steady_coordinates = self.nnormal_axes.coords_to_point(x, y)
            proj_tr = Dot(radius=0.05, color=RED).move_to(steady_coordinates)
            proj_trans_n.append(proj_tr)

        # get pdf and transformations of projections for simple regulation
        mean_s, var_s = mean(stable_states['simple']), variance(stable_states['simple'])
        pdf_s = self.snormal_axes.get_graph(lambda x_lvl: gauss(x_lvl, mean_s, var_s),
                                            color=BLUE,
                                            x_range=[-2, 8])
        proj_trans_s = []
        for x in stable_states['simple']:
            y = gauss(x, mean_s, var_s)
            steady_coordinates = self.snormal_axes.coords_to_point(x, y)
            proj_tr = Dot(radius=0.05, color=BLUE).move_to(steady_coordinates)
            proj_trans_s.append(proj_tr)

        # create scene
        self.play(FadeIn(self.title_dist))
        #self.add(self.nnormal_axes)
        #self.add(self.snormal_axes)

        # draw original projections
        # beta: [dot, lines, label]
        dots_proj = {'neg':VGroup(),'simple':VGroup()}
        lines_proj = {'neg':[],'simple':[]}
        for beta in self.beta_vary:
            self.play(Create(projections['neg'][beta][0]),Create(projections['simple'][beta][0]))
            self.play(Create(projections['neg'][beta][1]), Create(projections['simple'][beta][1]))
            dots_proj['neg'].add(projections['neg'][beta][0])
            dots_proj['simple'].add(projections['simple'][beta][0])
            lines_proj['neg'].append(projections['neg'][beta][1])
            lines_proj['simple'].append(projections['simple'][beta][1])
        # draw projection transformations
        self.play(Transform(dots_proj['neg'], VGroup(*proj_trans_n)))
        self.play(*[FadeOut(x) for x in lines_proj['neg']])
        self.play(Transform(dots_proj['simple'], VGroup(*proj_trans_s)),*[FadeOut(x) for x in lines_proj['simple']])
        self.play(Create(pdf_s),Create(pdf_n))

    def negative_auto_graph(self,beta=2,alpha=0.8,K=1,n=4):
        """
        plots production & degradation for negative autoregulation
        :return:
        """
        self.neg_protein = ProteinNegativeAuto(start=0,leaky_production_rate=0,
                                    max_production_rate=beta,degradation_rate=alpha)
        self.neg_protein.add_repressor(self.neg_protein,K,n)
        neg_graph = self.get_plots(self.neg_protein,axes=self.neg_axes)
        return neg_graph

    def simple_reg_graph(self,beta=2,alpha=0.8):
        """
        plots production & degradation for simple regulation
        :return:
        """
        self.simple_protein = ProteinSimpleReg(start=0,max_production_rate=beta,degradation_rate=alpha)
        simple_graph = self.get_plots(self.simple_protein,axes=self.simple_axes)
        return simple_graph

    def arrange_cells(self):
        cell1 = ImageMobject("assets/cell.png").scale(0.25).move_to(UR*2)
        cell2 = ImageMobject("assets/cell.png").scale(0.25).move_to(RIGHT*2)
        cell3 = ImageMobject("assets/cell.png").scale(0.25).move_to(DR*2)
        return Group(cell1,cell2,cell3)

    def simple_img(self):
        DNA = Rectangle(height=0.5, width=5, fill_opacity=0.1).set_color(BLACK).set_opacity(0.1).move_to(RIGHT)
        genes = VGroup(DNA)
        gene1 = Rectangle(height=0.5, width=1,fill_opacity=0.5).set_color(BLUE)
        gene2 = Rectangle(height=0.5, width=1,fill_opacity=0.5).set_color(PINK).next_to(gene1,buff=0.6)
        genes.add(gene1,gene2)
        protein = ImageMobject("assets/protein.png").scale(0.25).move_to(UP*1.5)
        text1, text2 = Text('gene X').move_to(gene1).set_color(BLACK).scale(0.5),\
                       Text('gene Y').move_to(gene2).set_color(BLACK).scale(0.5)
        prod_gene1 = Arrow(stroke_width=1,max_tip_length_to_length_ratio=0.1).set_color(BLACK)
        prod_gene1.put_start_and_end_on(gene1.get_top(),protein.get_bottom())
        activ_arrow = Arrow(protein.get_center(),gene2.get_top(),
                            stroke_width=1,color=BLACK,path_arc=-1,
                            max_tip_length_to_length_ratio=0.1)
        return Group(genes,text1,text2,protein,prod_gene1,activ_arrow)

    def negative_auto_img(self,negative=True):
        # production & degradation -> main functions of negative regulation
        gene1 = Rectangle(height=0.5, width=2,fill_opacity=0.5).set_color(BLUE)
        protein = ImageMobject("assets/protein.png").scale(0.25).move_to(UP*1.5)
        text1 = Text('gene X').move_to(gene1).set_color(BLACK).scale(0.5)
        prod_gene1 = Arrow(stroke_width=1,max_tip_length_to_length_ratio=0.1).set_color(BLACK)
        prod_gene1.put_start_and_end_on(gene1.get_top(),protein.get_center())
        if negative:
            back_arrow = Arrow(protein.get_center(), gene1.get_center(),
                               stroke_width=1, color=BLACK, path_arc=-1,
                               max_tip_length_to_length_ratio=0.1,tip_shape=InhibitionArrowTip)
        else:
            back_arrow = Arrow(protein.get_center(), gene1.get_center(),
                               stroke_width=1, color=BLACK, path_arc=-1,
                               max_tip_length_to_length_ratio=0.1)
        return Group(gene1,text1,protein,prod_gene1,back_arrow)

    def get_projection(self,protein=None,axes=None,label='Y_{st}'):
        if not axes:
            if type(protein) == ProteinNegativeAuto:
                axes = self.neg_axes
            elif type(protein) == ProteinSimpleReg:
                axes = self.simple_axes
        Xst = bisection(lambda x_lvl: protein.production_rate(x_lvl),
                        lambda x_lvl: protein.degradation_rate(x_lvl),
                        0, 10)
        Yst_prod = protein.production_rate(Xst)
        steady_coordinates = axes.coords_to_point(Xst, Yst_prod)
        project_coordinates = axes.coords_to_point(Xst, 0)
        dot = Dot(radius=0.05, color=BLACK).move_to(steady_coordinates)
        lines = DashedVMobject(Line(steady_coordinates, project_coordinates, color=BLACK))
        Xst_label = MathTex(label, color=BLACK).next_to(project_coordinates + DOWN).scale(0.8)
        return [dot, lines, Xst_label]

    def get_plots(self,protein,axes=None,plot_types=set(['prod','deg','max','steady','K']),col = BLUE,xr=[0,7]):
        res = {}
        if not axes:
            axes = self.axes
        if 'prod' in plot_types:
            prod_f = lambda x: protein.production_rate(x)
            prod_graph = axes.get_graph(lambda x_lvl: prod_f(x_lvl),
                                             color=col,
                                             x_range=xr)
            prod_label = axes.get_graph_label(
                prod_graph, Tex(r"production"), x_val=xr[1], direction= UP /3
            ).scale(0.5)
            res['prod'] = [prod_graph, prod_label]
        if 'deg' in plot_types:
            deg_f = lambda x: protein.degradation_rate(x)
            deg_graph = axes.get_graph(lambda x_lvl: deg_f(x_lvl),
                                            color=RED,
                                            x_range=xr)
            deg_label = axes.get_graph_label(
                deg_graph, Tex(r"degradation"), x_val=xr[1], direction=UP / 3
            ).scale(0.5)
            res['deg'] = [deg_graph, deg_label]
        if 'max' in plot_types:
            max_plot = axes.get_graph(lambda x_lvl: protein.production_rates['max'],
                                       color=col,
                                       x_range=xr)
            max_graph = DashedVMobject(max_plot)
            max_label = axes.get_graph_label(
                max_plot, Tex(r"max prod."), x_val=-3, direction= RIGHT
            ).scale(0.5)
            res['max']=[max_graph,max_label]
        if 'K' in plot_types:
            # K is a binding property of promoter
            if type(protein) == ProteinSimpleReg:
                res['K'] = [Dot().scale(0.01),Dot().scale(0.01)]
            else:
                xK = protein.r_kinetics['K']
                yK = protein.production_rate(xK)
                K_graph = axes.get_graph(lambda x_lvl: xK,
                                              color=GREEN,
                                              x_range=[0, xK])
                K_label = axes.get_graph_label(
                    K_graph, Tex(r"K"), x_val=xK, direction=UP / 2
                )
                res['K'] = [K_graph,K_label]
        if 'steady' in plot_types:
            res['steady'] = self.get_projection(protein,axes)
        return res

    def brace_label(self,objects,label, brace_direction=RIGHT):
        """

        :param brace_direction: UP / DOWN / LEFT / RIGHT
        :param label: str, how to label group of objects
        :param objects: list of objects to label
        :return: BraceText
        """
        kwargs = {"buff":0}
        group = VGroup(*objects)
        b = Brace(group,direction=brace_direction,color=BLACK)
        t = b.get_text(label,**kwargs).set_color(BLACK).scale(0.5)
        return VGroup(b,t)


    def Appear(self,obj,lag = 0.5):
        """

        :param obj: Mobject to draw and remove
        :return: None, plays the object instead
        """
        self.play(FadeIn(obj))
        self.wait(lag)
        self.play(FadeOut(obj))


def bisection(f,g, a, b):
  # returns x for which f and g functions intersect
    eps = 10 ** -10
    h = lambda x: f(x)-g(x)
    ha = h(a)
    hb = h(b)
    if ha * hb > 0:
        raise ValueError("Bad input")

    for i in range(1000):
        ha = h(a)
        midpoint = (a + b) / 2
        hm = h(midpoint)

        if abs(hm) < eps:
            return midpoint

        if hm * ha >= 0:
            a = midpoint
        else:
            b = midpoint

    raise RuntimeError("Out of iterations")

class InhibitionArrowTip(ArrowTip, Rectangle):
    def __init__(self, length=0.35, **kwargs):
        Rectangle.__init__(self, height=0,width=0.00000001, **kwargs)
        self.width = length
        self.stretch_to_fit_height(length)