from manim import *
from bio import *

config.background_color = WHITE

class HillExplanation(Scene):

    def construct(self):
        self.initialize_proteins()
        self.initialize_scene()
        self.X_graph_func()
        self.wait(1)
        # self.max_lvl()
        # self.wait(1)

    def initialize_proteins(self,X_beta=0.1, Y_beta = 1,K=0.5,n=4):
        # assumes there is no degradation and no additional regulation for X
        self.X = Protein(beta=X_beta)
        self.Y = Protein(beta=Y_beta)
        self.K = K # activation constant
        self.n = n # steepness of hill function

    def max_lvl(self):
        # any protein has maximum level of production
        X_max = self.X_axes.get_graph(lambda t: self.X.beta, x_range=[0,10], color=RED)
        Xmax_label = self.X_axes.get_graph_label(
            X_max, MathTex(r"X_{max}"), x_val=10
        ).scale(0.7)
        Xmax_line = VGroup(X_max,Xmax_label)

        beta = self.dYdt_axes.get_graph(lambda x_lvl: self.Y.beta, x_range=[0,3], color=TEAL_E)
        beta_label = self.dYdt_axes.get_graph_label(
            beta, MathTex(r"\beta"), x_val=3
        ).scale(0.7)
        beta_line = VGroup(beta,beta_label)
        explanation_1 = MathTex(r"\text{Maximum \; level \,}",
                                r"\text{of \, X}",
                                r"\\ \text{that \, can \, be \, produced}",
                                color=BLACK).scale(0.5).next_to(Xmax_label, RIGHT)
        explanation_2 = MathTex(r"\text{Maximum \; rate} \\",
                                r"\text{of \,production}",
                                r"\text{ of \, Y}",
                                color=BLACK).scale(0.5).next_to(beta_label, RIGHT)

        self.play(Create(Xmax_line),Create(beta_line),Create(explanation_1),Create(explanation_2))
        self.play(Indicate(explanation_1[0],color=RED),
                           Indicate(explanation_2[0],color=TEAL_E))

    def initialize_scene(self):
        # what the scene is about
        name = Tex("Why Hill function?", color=BLACK).scale(2).move_to(UP)
        reg_text = Tex("X→Y", color=BLACK).move_to(DOWN)
        self.play(Create(name), Create(reg_text))
        self.wait(1)
        self.play(FadeOut(name),reg_text.animate.move_to(UP*3 + LEFT * 5))

        # number plane, graphs

        self.number_plane = NumberPlane(
            x_range=[-10, 10, 0.5],
            y_range=[-10, 10, 0.5],
            background_line_style={
                "stroke_color": BLACK,
                "stroke_width": 1,
                "stroke_opacity": 0.4
            })
        self.number_plane.coords_to_point(0)
        self.dYdt_graph_scene()
        self.X_vs_t_graph_scene()
        self.graph_group = VGroup(self.dYdt_graph,self.X_graph)
        self.graph_group.arrange_in_grid(2, 1).scale(0.5)
        self.X_axes.coords_to_point(0)
        self.dYdt_axes.coords_to_point(0)
        self.play(Create(self.graph_group))
        #self.play(Create(self.number_plane))

    def dYdt_graph_scene(self):
        self.dYdt_axes = Axes(
            x_range=[0, 3, 0.5],
            y_range=[0, 1.5, 0.5],
            axis_config={"color": BLACK, "include_tip": False},

        )
        self.dYdt_axes.coords_to_point(0)
        y_label = self.dYdt_axes.get_y_axis_label(MathTex(r"\frac{dY}{dt}"),
                                                  edge=LEFT, direction=LEFT, buff=0.4).set_color(BLACK)
        x_label = self.dYdt_axes.get_x_axis_label("X").set_color(BLACK)
        self.dYdt_labels = VGroup(x_label, y_label)
        self.dYdt_graph = VGroup(self.dYdt_axes,self.dYdt_labels)


    def X_vs_t_graph_scene(self):
        self.X_axes = Axes(
            x_range=[0, 10, DELTA_T],
            y_range=[0, 3, 0.5],
            axis_config={"color": BLACK, "include_tip": False},
        )
        self.X_axes.coords_to_point(0)
        y_label = self.X_axes.get_y_axis_label(MathTex(r"X \\ \text{(activator of Y)}"),
                                               edge=LEFT, direction=LEFT, buff=0.4).set_color(BLACK).move_to(LEFT)
        x_label = self.X_axes.get_x_axis_label("t").set_color(BLACK)
        self.X_labels = VGroup(x_label, y_label)
        self.X_graph = VGroup(self.X_axes, self.X_labels)

    def dYdt_graph_func(self):
        pass

    def X_graph_func(self,t=10):
        dot_collection = VGroup()
        for time, val in self.X.level_record.items():
            dot = Dot(color=RED).move_to(self.X_axes.coords_to_point(time, val))
            self.X_axes.add(dot)
            dot_collection.add(dot)

    def X_lvls_updater(self):
        """
        self.setup_axes()
        graph = self.get_graph(lambda x: 0.1 * (x + 7) * (x - 2) * (x - 7), x_min=-15, x_max=15)

        start_x = -10
        end_x = 10
        tracker = ValueTracker(start_x)  # starting point of x

        ctp = self.coords_to_point
        dot_x = tracker.get_value
        func = graph.underlying_function
        moving_dot = always_redraw(lambda: Dot(ctp(dot_x(), func(dot_x())), radius=0.1, color=RED, ))

        self.add(graph, moving_dot)
        self.play(tracker.set_value, end_x, run_time=5)
        self.wait()

        :return:
        """
