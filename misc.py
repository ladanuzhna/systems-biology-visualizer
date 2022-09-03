from manim import *

def simple_regulation(self):
    # returns X->Y Mobject
    v1 = Dot(radius=0.5).set_color(RED_A)
    v2 = Dot(radius=0.5).set_color(RED_E).to_corner(v1.get_center() + [10 ** (-50), 10 ** (-50), 0])
    edge = Arrow(max_stroke_width_to_length_ratio=0.5, path_arc=2, color=RED)
    edge.put_start_and_end_on(v1.get_center(), v2.get_start())
    self.play(FadeIn(v1), Create(edge))
    self.play(Create(v2))
    return VGroup(v1, v2, edge)

class BraceExample(Scene):

    def construct(self):
        t1 = Text("text 1")
        t2 = Text("text 2")
        group = VGroup(t1,t2).arrange(DOWN)
        b = BraceText(group,"group of text",brace_direction=RIGHT)
        self.add(group)
        self.add(b)