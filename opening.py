from manim import *
import networkx as nx

config.background_color = WHITE

class OpeningManim(Scene):
    def construct(self):
        self.construct_intro()
        self.construct_graph()

    def construct_intro(self):
        title = Tex(r"Simple regulation X->Y ",color=BLACK)
        basel = MathTex(r"\frac{dY}{dt} = \beta - \alpha*Y",color=BLACK)
        VGroup(title, basel).arrange(DOWN)
        self.play(
            Write(title),
            FadeIn(basel, shift=DOWN),run_time=1
        )

        transform_title = Tex("X->Y",color=BLACK)
        transform_title.to_corner(UP + LEFT)
        self.play(
            Transform(title, transform_title),
            LaggedStart(*[FadeOut(obj, shift=DOWN) for obj in basel]),
        )

    def construct_graph(self, n_nodes=14, connectivity = 1):
        nxgraph = nx.erdos_renyi_graph(n_nodes, connectivity)
        G = Graph.from_networkx(nxgraph,
                                layout="spring",layout_scale=4).set_color(BLACK)
        self.play(Create(G))
        unpack_graph = []
        for ind, v in enumerate(G.vertices):
            unpack_graph.append(G[v].animate.move_to(5 * RIGHT * np.cos(ind / 7 * PI) +
                                         3 * UP * np.sin(ind / 7 * PI)))
            unpack_graph.append(Transform((G.vertices[ind],G.vertices[ind].copy().scale(1.5))))
        unpack_graph.append(G.vertices[0].animate.set_color(RED))
        unpack_graph.append(G.vertices[5].animate.set_color(RED))
        unpack_graph.append(G.edges[(0,5)].animate.set_color(RED))
        self.play(*unpack_graph)
        #moved_G = Graph([0,5,4],[(5,4)]).set_color(RED).to_corner(LEFT)
        self.play(*[G._remove_edges_animation(*[e for e in G.edges if e != (0,5)]),
                    G._remove_vertices_animation(*[v for v in range(n_nodes) if v != 0 and v != 5])])
        self.play(Transform(G,self.simple_regulation()))
        self.wait(0.3)

    def simple_regulation(self,color=RED):
        # returns X->Y Mobject
        v1 = Dot(radius=0.2).set_color(color)
        v2 = Dot(radius=0.2).set_color(color).to_corner(v1.get_center()+[-0.05,0.001,0])
        edge = Arrow(max_stroke_width_to_length_ratio=0.5,
                     max_tip_length_to_length_ratio=0.3,
                     color=color)
        edge.put_start_and_end_on(v1.get_center(), v2.get_center())
        return VGroup(v1,v2,edge).scale(0.5)

    def positive_double_regulation(self):
        # returns X<->Y Mobject
        pass

    def negative_double_regulation(self):
        # returns X|-|Y Mobject
        pass

    def negative_self_regulation(self):
        # returns X|-|X Mobject
        pass

    def positive_self_regulation(self):
        # returns X<->X Mobject
        pass