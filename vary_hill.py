from manim import *
from bio import *

config.background_color = WHITE

# graph that shows varying parameters of Hill function

# TODO: show upregulation X->Y graph
# TODO: do the same simulation for downregulation

class HillUpregulation(MovingCameraScene):

    def construct(self,x_graphing_range=[0,1]):
        self.range = x_graphing_range
        self.graph_scene()
        self.beta_var,self.K_var,self.n_var =1,0.5,4
        self.on_screen_var_beta = Variable(self.beta_var,
                                      MathTex(r"\beta"),
                                      num_decimal_places=1).set_color(BLACK).move_to(LEFT*2 + DOWN*3.5).scale(0.7)
        self.on_screen_var_K = Variable(self.K_var,
                                      MathTex(r"K"),
                                      num_decimal_places=1).set_color(BLACK).move_to( DOWN*3.5).scale(0.7)
        self.on_screen_var_n = Variable(self.n_var,
                                      MathTex(r"n"),
                                      num_decimal_places=0).set_color(BLACK).move_to(RIGHT*2 + DOWN*3.5).scale(0.7)
        self.eqn, self.focus = self.simple_reg_eqn()
        reg_text = Tex("X->Y", color=BLACK).scale(2)
        self.play(Create(reg_text,lag_ratio=0.1),Create(self.number_plane))
        self.wait(1)
        self.play(Create(self.eqn),reg_text.animate.scale(0.4).move_to(UP*3 + LEFT*4))
        self.beta_graph()
        self.K_graph()
        self.n_graph()
        self.play(FadeOut(self.on_screen_var_beta),
                  FadeOut(self.on_screen_var_K),
                  FadeOut(self.on_screen_var_n))

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
        self.axes = Axes(
            x_range=[0, 1.2, 0.5],
            y_range=[0, 1, 0.5],
            axis_config={"color": BLACK, "include_tip": False},
        )
        self.axes.coords_to_point(0)

        y_label = self.axes.get_y_axis_label(MathTex(r"\frac{dY}{dt}"), edge=LEFT, direction=LEFT, buff=0.4).set_color(
            BLACK)
        x_label = self.axes.get_x_axis_label("X").set_color(BLACK)
        self.labels = VGroup(x_label, y_label)

    def beta_graph(self):
        vary = lambda X, beta, K, n: Protein(beta=beta).levels_vs_activator(X, K, n)
        # Create Graph
        graph = self.axes.get_graph(lambda x_lvl: vary(x_lvl, beta=self.beta_var, K=self.K_var, n=self.n_var), color=MAROON,
                               x_range=self.range)
        graph_label = self.axes.get_graph_label(
            graph, MathTex(r"\beta="+str(self.beta_var)), x_val=1.05
        ).scale(0.7)
        self.play(Create(self.focus["beta"]),DrawBorderThenFill(self.axes),Write(graph_label),
                  Write(self.labels),Create(graph),
                  Write(self.on_screen_var_K),Write(self.on_screen_var_beta),Write(self.on_screen_var_n))
        graphs = [graph]
        graph_labels = [graph_label]
        colors = [ MAROON_E, MAROON_B]
        for i,beta in enumerate([0.5,0.1]):
            if i >= len(colors):
                raise Exception("Not enough colors for the graphs generated")
            graph_ = self.axes.get_graph(lambda x_lvl: vary(X=x_lvl,beta=beta, K=self.K_var, n=self.n_var),
                                       color=colors[i],x_range = self.range )
            graph_label = self.axes.get_graph_label(
                graph_, MathTex(r"\beta=}" + str(beta)), x_val=1.05
            ).scale(0.7)
            self.beta_var = beta
            copy = graph.copy()
            self.play(Transform(copy, graph_),Create(graph_label),
                self.on_screen_var_beta.tracker.animate.set_value(beta))
            graph_labels.append(graph_label)
            graphs.append(graph_)
            graphs.append(copy)
            self.wait(0.4)

        group = VGroup(*graphs[1:],*graph_labels)
        self.beta_var =  1
        self.play(FadeOut(group),graphs[0].animate.set_color(BLUE),
                  self.on_screen_var_beta.tracker.animate.set_value(self.beta_var))

    def K_graph(self):
        vary = lambda X, beta, K, n: Protein(beta=beta).levels_vs_activator(X, K, n)
        # Create Graph
        graph = self.axes.get_graph(lambda x_lvl: vary(x_lvl, beta=self.beta_var, K=self.K_var, n=self.n_var),
                                    color=BLUE,
                                    x_range=self.range)
        graph_label = self.axes.get_graph_label(
            graph, MathTex(r"K=" + str(self.K_var)), x_val=1.05, direction=UP / 2
        ).scale(0.7)
        self.play(
            ReplacementTransform(self.focus["beta"], self.focus["K"])
        )
        self.play(Write(graph_label),Create(graph))
        graphs = [graph]
        graph_labels = [graph_label]
        colors = [TEAL_E, GREEN_C]
        for i, K in enumerate([1, 0.1]):
            if i >= len(colors):
                raise Exception("Not enough colors for the graphs generated")
            graph_ = self.axes.get_graph(lambda x_lvl: vary(X=x_lvl, beta=self.beta_var, K=K, n=self.n_var),
                                         color=colors[i], x_range=self.range)
            graph_label = self.axes.get_graph_label(
                graph_, MathTex(r"K=" + str(K)), x_val=1.05, direction=UP / 2
            ).scale(0.7)
            self.K_var = K
            copy = graph.copy()
            self.play(Transform(copy, graph_), Create(graph_label),
                      self.on_screen_var_K.tracker.animate.set_value(K))
            graph_labels.append(graph_label)
            graphs.append(graph_)
            graphs.append(copy)
            self.wait(0.4)

        group = VGroup(*graphs[1:], *graph_labels)
        self.K_var = 0.5
        self.play(FadeOut(group), graphs[0].animate.set_color(GREEN),
                  self.on_screen_var_K.tracker.animate.set_value(self.K_var))

    def n_graph(self):
        vary = lambda X, beta, K, n: Protein(beta=beta).levels_vs_activator(X, K, n)
        # Create Graph
        graph = self.axes.get_graph(lambda x_lvl: vary(x_lvl, beta=self.beta_var, K=self.K_var, n=self.n_var),
                                    color=GREEN,
                                    x_range=self.range)
        graph_label = self.axes.get_graph_label(
            graph, MathTex(r"n=" + str(self.n_var)), x_val=1.05, direction=UP / 2
        ).scale(0.7)
        self.play(
            ReplacementTransform(self.focus["K"], self.focus["n"])
        )
        self.play(Write(graph_label), Create(graph))
        graphs = [graph]
        graph_labels = [graph_label]
        colors = [GREEN_E, TEAL_D,YELLOW_E,BLUE_E ]
        for i, n in enumerate([3,2,1]):
            if i >= len(colors):
                raise Exception("Not enough colors for the graphs generated")
            graph_ = self.axes.get_graph(lambda x_lvl: vary(X=x_lvl, beta=self.beta_var, K=self.K_var, n=n),
                                         color=colors[i], x_range=self.range)
            graph_label = self.axes.get_graph_label(
                graph_, MathTex(r"n=" + str(n)), x_val=1.05, direction=UP / 2
            ).scale(0.7)
            self.n_var = n
            copy = graph.copy()
            self.play(Transform(copy, graph_), Create(graph_label),
                      self.on_screen_var_n.tracker.animate.set_value(n))
            graph_labels.append(graph_label)
            graphs.append(graph_)
            graphs.append(copy)
            self.wait(0.4)

        group = VGroup(*graphs[1:], *graph_labels)
        self.n_var = 4
        self.play(FadeOut(group), FadeOut(self.focus["n"]),graphs[0].animate.set_color(GREEN),
                  self.on_screen_var_n.tracker.animate.set_value(self.n_var))

    def simple_reg_eqn(self):
        # TODO: make sure it is downregulation / upregulation
        eqn = MathTex(r"\frac{dY}{dt}=",r"\beta",
                     r"{[X]",r"\raise0.5ex\hbox{n}",r"\over",
                      r"[K]",r"\raise0.5ex\hbox{n}", "+", r"[X]",r"\raise0.5ex\hbox{n}}").set_color(BLACK).scale(0.7).move_to(UP*3+LEFT)
        focus_on = {"K":SurroundingRectangle(eqn[5], buff = .1).set_color(BLUE),
                    "n":VGroup(SurroundingRectangle(eqn[3], buff=.1).set_color(GREEN),
                         SurroundingRectangle(eqn[6], buff=.1).set_color(GREEN),
                         SurroundingRectangle(eqn[9], buff=.1).set_color(GREEN)),
                    "beta":SurroundingRectangle(eqn[1], buff=.1).set_color(MAROON)}
        return eqn,focus_on




"""
circ = Circle().scale(1.5)
        circ_ref = circ.copy()
        circ.apply_complex_function(
            lambda x: np.exp(x*1j)
        )
        t = ValueTracker(0)
        circ.add_updater(
            lambda x: x.become(circ_ref.copy().apply_complex_function(
                lambda x: np.exp(x+t.get_value()*1j)
            )).set_color(BLUE)
        )
        self.add(circ_ref)
        self.play(TransformFromCopy(circ_ref, circ))
        self.play(t.animate.set_value(TAU), run_time=3)
"""
