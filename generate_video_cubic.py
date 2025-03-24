from gtts import gTTS
from sympy import symbols, solve, sympify, diff, discriminant
import numpy as np
from manim import *
from pydub import AudioSegment

def generate_video_cubic(equation):
    # Parse and analyze the cubic equation
    x = symbols('x')
    eq = sympify(equation)
    
    # Extract coefficients (a*x^3 + b*x^2 + c*x + d)
    coeffs = eq.as_poly(x).all_coeffs()
    while len(coeffs) < 4:  # Pad with zeros if lower degree
        coeffs.insert(0, 0)
    a, b, c, d = coeffs
    
    # Handle degenerate case (not a cubic)
    is_cubic = (a != 0)
    
    # Find roots
    roots = list(solve(eq, x))
    real_roots = [float(root.evalf()) for root in roots if root.is_real]
    complex_roots = [root for root in roots if not root.is_real]
    from collections import Counter
    roots_dict = dict(Counter(roots))
    
    # Find derivative for critical points
    derivative = diff(eq, x)
    critical_points = solve(derivative, x)
    real_critical_points = [float(cp.evalf()) for cp in critical_points if cp.is_real]
    
    # Find second derivative for inflection point
    second_derivative = diff(derivative, x)
    inflection_points = solve(second_derivative, x)
    real_inflection_points = [float(ip.evalf()) for ip in inflection_points if ip.is_real]
    
    # Dynamic axis scaling
    key_points = real_roots + real_critical_points + real_inflection_points
    if not key_points:
        x_range = [-5, 5]
        y_range = [-8, 8]
    else:
        x_min, x_max = min(key_points) - 2, max(key_points) + 2
        x_range = [x_min, x_max]
        y_vals = [float(eq.subs(x, pt).evalf()) for pt in np.linspace(x_min, x_max, 100)]
        y_range = [min(y_vals) - 2, max(y_vals) + 2]
    
    # Create verbal equation for narration
    wordEquation = ''
    if a == 1:
        wordEquation += "x cubed "
    elif a == -1:
        wordEquation += "minus x cubed "
    elif a < 0:
        wordEquation += f"minus {abs(a)} x cubed "
    elif a > 0:
        wordEquation += f"{a} x cubed "
    
    if b == 1:
        wordEquation += "plus x squared "
    elif b == -1:
        wordEquation += "minus x squared "
    elif b > 0:
        wordEquation += f"plus {b} x squared "
    elif b < 0:
        wordEquation += f"minus {abs(b)} x squared "
    
    if c == 1:
        wordEquation += "plus x "
    elif c == -1:
        wordEquation += "minus x "
    elif c > 0:
        wordEquation += f"plus {c} x "
    elif c < 0:
        wordEquation += f"minus {abs(c)} x "
    
    if d > 0:
        wordEquation += f"plus {d}"
    elif d < 0:
        wordEquation += f"minus {abs(d)}"
    elif len(coeffs) == 1 and a == 0:
        wordEquation = "zero"
    
    # Generate script parts
    intro = f"Welcome to Animation. I am your instructor Jojo. You have given the equation {wordEquation}."
    
    # Format roots for narration
    roots_str = ""
    if not is_cubic:
        roots_str = "This is not a cubic equation but a lower-degree polynomial."
    elif len(real_roots) == 3 and len(roots_dict) == 1:  # Triple root
        roots_str = f"It has a triple real root at x equals {float(real_roots[0]):.2f}."
    elif len(real_roots) == 3:  # Three distinct roots
        roots_str = f"It has three real roots at x equals {float(real_roots[0]):.2f}, {float(real_roots[1]):.2f}, and {float(real_roots[2]):.2f}."
    elif len(real_roots) == 2:  # Double root + distinct root
        repeated_root = [r.evalf() for r in roots if roots_dict[r] == 2][0]
        distinct_root = [r.evalf() for r in roots if roots_dict[r] == 1][0]
        roots_str = f"It has a double root at x equals {float(repeated_root):.2f} and a distinct root at x equals {float(distinct_root):.2f}."
    elif len(real_roots) == 1:
        roots_str = f"It has one real root at x equals {float(real_roots[0]):.2f} and two complex roots."
    else:
        roots_str = f"It has no real roots, only three complex roots."
    
    # Differentiation functions
    def cubic_function(x):
        return float((a * x**3 + b * x**2 + c * x + d).evalf())
    
    def diff_function(x):
        return float((3 * a * x**2 + 2 * b * x + c).evalf())
    
    def second_diff_function(x):
        return float((6 * a * x + 2 * b).evalf())
    
    # Format critical points for narration
    critical_str = ""
    if not is_cubic or a == 0 and b == 0:  # Linear or constant
        critical_str = "There are no turning points for this non-cubic equation."
    elif len(real_critical_points) == 2:
        cp1_x, cp2_x = sorted(real_critical_points)
        cp1_y = cubic_function(cp1_x)
        cp2_y = cubic_function(cp2_x)
        if second_diff_function(cp1_x) > 0:
            critical_str = f"There are two turning points: a local minimum at x equals {cp1_x:.2f} with f(x) equals {cp1_y:.2f}, and a local maximum at x equals {cp2_x:.2f} with f(x) equals {cp2_y:.2f}."
        else:
            critical_str = f"There are two turning points: a local maximum at x equals {cp1_x:.2f} with f(x) equals {cp1_y:.2f}, and a local minimum at x equals {cp2_x:.2f} with f(x) equals {cp2_y:.2f}."
    elif len(real_critical_points) == 1:
        cp_x = real_critical_points[0]
        cp_y = cubic_function(cp_x)
        critical_str = f"There is one critical point at x equals {cp_x:.2f} with f(x) equals {cp_y:.2f}, a horizontal inflection due to a repeated root."
    else:
        critical_str = "There are no real turning points."
    
    # Format inflection point for narration
    inflection_str = ""
    if not is_cubic or a == 0 and b == 0:
        inflection_str = "There is no inflection point for this non-cubic equation."
    elif len(real_inflection_points) == 1:
        inf_x = real_inflection_points[0]
        inf_y = cubic_function(inf_x)
        inflection_str = f"There is an inflection point at x equals {inf_x:.2f} with f(x) equals {inf_y:.2f}, where the concavity changes."
    else:
        inflection_str = "Unexpectedly, no real inflection points were found."
    
    # Describe end behavior
    if a > 0:
        behavior = "As x approaches positive infinity, f(x) approaches positive infinity. As x approaches negative infinity, f(x) approaches negative infinity."
    elif a < 0:
        behavior = "As x approaches positive infinity, f(x) approaches negative infinity. As x approaches negative infinity, f(x) approaches positive infinity."
    elif b > 0:
        behavior = "This is a quadratic with end behavior rising to positive infinity on both sides."
    elif b < 0:
        behavior = "This is a quadratic with end behavior falling to negative infinity on both sides."
    else:
        behavior = "This is a linear or constant function with no cubic behavior."
    
    # Final narration for drawing
    draw = "Let's draw and examine this equation graph!"
    slope = "Now visualize the slope of the graph."
    script = [intro, draw, behavior, roots_str, critical_str, inflection_str, slope]
    
    # Generate audio files
    audio_paths = []
    for i, text in enumerate(script):
        audio_path = f"audio_{i+1}.mp3"
        audio = gTTS(text=text, lang='en')
        audio.save(audio_path)
        audio_paths.append(audio_path)
    
    # Load audio segments for timing
    audio_segments = [AudioSegment.from_mp3(path) for path in audio_paths]
    audio_lengths = [len(segment)/1000 for segment in audio_segments]  # Convert to seconds

    class CubicFunctionExplanation(MovingCameraScene):
        def construct(self):
            # Configuration
            axes_config = {
                "x_range": [x_range[0], x_range[1], 1],
                "y_range": [y_range[0], y_range[1], max(2, (y_range[1] - y_range[0]) / 4)],
                "axis_config": {"color": BLUE},
                "x_axis_config": {"numbers_to_include": np.arange(int(x_range[0]), int(x_range[1]) + 1, 2), "font_size": 18},
                "y_axis_config": {"numbers_to_include": np.arange(int(y_range[0]), int(y_range[1]) + 1, max(2, int((y_range[1] - y_range[0]) / 4))), "font_size": 18},
            }
            
            # Create coordinate system
            axes = Axes(**axes_config)
            x_label = MathTex("x", font_size=20).next_to(axes.x_axis.get_end(), RIGHT)
            y_label = MathTex("f(x)", font_size=20).next_to(axes.y_axis.get_end(), UP)
            axes_labels = VGroup(x_label, y_label)
            
            # Title
            title = Text("Understanding Cubic Equations" if is_cubic else "Understanding Polynomials", font_size=48)
            subtitle = Text(f"f(x) = {equation}", font_size=36)
            subtitle.next_to(title, DOWN)
            self.add_sound(audio_paths[0])
            self.play(AddTextLetterByLetter(title), run_time=audio_lengths[0]*0.6)
            self.play(AddTextLetterByLetter(subtitle), run_time=audio_lengths[0]*0.4)
            self.wait(0.5)
            self.play(title.animate.scale(1.2), subtitle.animate.scale(1.2), run_time=0.5)
            self.play(title.animate.to_edge(UP).scale(0.8333).set_opacity(0),
                      subtitle.animate.to_edge(DOWN).scale(0.8333).set_opacity(0), run_time=1)
            
            # Display coordinate system
            self.add_sound(audio_paths[1])
            self.play(Create(axes), Write(axes_labels))
            self.wait(1)
            
            # Create function graph
            function_label = MathTex(f"f(x) = {equation}")
            function_label.to_corner(UL)
            graph = axes.plot(cubic_function, color=GREEN)
            self.play(Create(graph))
            self.play(Write(function_label))
            self.wait(1)
            self.play(FadeOut(function_label))
            
            # Create arrows showing behavior at infinities
            if a > 0:
                right_arrow = Arrow(axes.c2p(x_range[1] - 1, y_range[1] - 1), axes.c2p(x_range[1], y_range[1]), color=YELLOW)
                left_arrow = Arrow(axes.c2p(x_range[0] + 1, y_range[0] + 1), axes.c2p(x_range[0], y_range[0]), color=YELLOW)
                right_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
                left_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
            elif a < 0:
                right_arrow = Arrow(axes.c2p(x_range[1] - 1, y_range[0] + 1), axes.c2p(x_range[1], y_range[0]), color=YELLOW)
                left_arrow = Arrow(axes.c2p(x_range[0] + 1, y_range[1] - 1), axes.c2p(x_range[0], y_range[1]), color=YELLOW)
                right_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
                left_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
            else:  # Quadratic or lower
                if b > 0:
                    right_arrow = Arrow(axes.c2p(x_range[1] - 1, y_range[1] - 1), axes.c2p(x_range[1], y_range[1]), color=YELLOW)
                    left_arrow = Arrow(axes.c2p(x_range[0] + 1, y_range[1] - 1), axes.c2p(x_range[0], y_range[1]), color=YELLOW)
                    right_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
                    left_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
                elif b < 0:
                    right_arrow = Arrow(axes.c2p(x_range[1] - 1, y_range[0] + 1), axes.c2p(x_range[1], y_range[0]), color=YELLOW)
                    left_arrow = Arrow(axes.c2p(x_range[0] + 1, y_range[0] + 1), axes.c2p(x_range[0], y_range[0]), color=YELLOW)
                    right_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
                    left_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
                else:
                    right_arrow = left_arrow = right_text = left_text = VGroup()  # No arrows for linear/constant
            
            self.add_sound(audio_paths[2])
            self.play(Create(right_arrow), Create(left_arrow), Write(right_text), Write(left_text), run_time=audio_lengths[2])
            self.play(FadeOut(right_arrow), FadeOut(left_arrow), FadeOut(right_text), FadeOut(left_text))
            
            # Roots
            roots_title = Text("Roots", font_size=32).to_corner(UL)
            self.play(Write(roots_title))
            self.add_sound(audio_paths[3])
            root_dots = []
            root_labels = []
            graph_group = VGroup(axes, x_label, y_label, graph)
            if real_roots:
                for i, root in enumerate(sorted(set(real_roots))):
                    dot = Dot(axes.c2p(root, 0), color=RED)
                    multiplicity = roots_dict.get(root, 1)
                    label = MathTex(f"x_{i+1} = {root:.2f}" + (f" (m={multiplicity})" if multiplicity > 1 else ""), color=RED, font_size=24).next_to(dot, DOWN)
                    root_dots.append(dot)
                    root_labels.append(label)
                    self.play(self.camera.frame.animate.set_width(2).move_to(dot.get_center()), Create(dot), Write(label), run_time=1)
                    self.wait()
                    self.play(self.camera.frame.animate.set_width(graph_group.width * 1.3).move_to(axes.get_center()), run_time=0.5)
            else:
                no_roots_text = Text("No real roots for this equation", color=RED, font_size=24).next_to(axes.c2p(0, 0), DOWN)
                self.play(Write(no_roots_text))
                self.wait(1)
                self.play(FadeOut(no_roots_text))
            self.play(*[FadeOut(obj) for obj in root_dots + root_labels + [roots_title]])
            
            # Turning Points
            turning_points_title = Text("Turning Points", font_size=32).to_corner(UL)
            self.add_sound(audio_paths[4])
            self.play(Write(turning_points_title))
            if real_critical_points:
                for tp in sorted(real_critical_points):
                    y_tp = cubic_function(tp)
                    tp_dot = Dot(axes.c2p(tp, y_tp), color=YELLOW)
                    tp_label = MathTex(f"({tp:.2f}, {y_tp:.2f})", color=YELLOW, font_size=24).next_to(tp_dot, UP)
                    tangent_line = axes.plot(lambda x: y_tp, x_range=[tp - 1, tp + 1], color=YELLOW)
                    self.play(Create(tp_dot), Write(tp_label))
                    self.wait()
                    self.play(Create(tangent_line))
                    self.wait(0.3)
                    sdd = second_diff_function(tp)
                    if sdd > 0:
                        max_min_label = Text("Local Minimum", font_size=20, color=YELLOW)
                    elif sdd < 0:
                        max_min_label = Text("Local Maximum", font_size=20, color=YELLOW)
                    else:
                        max_min_label = Text("Horizontal Inflection", font_size=20, color=YELLOW)
                    max_min_label.next_to(tp_label, UP)
                    self.play(Write(max_min_label))
                    self.wait(1)
                    self.play(FadeOut(tangent_line), FadeOut(tp_dot), FadeOut(tp_label), FadeOut(max_min_label))
            else:
                turning_point_label = Text('No turning points for this equation.').to_edge(DOWN)
                self.play(Write(turning_point_label), run_time=audio_lengths[4])
                self.wait(1)
                self.play(FadeOut(turning_point_label))
            self.play(FadeOut(turning_points_title))
            
            # Inflection Point
            inflection_title = Text("Inflection Point", font_size=32).to_corner(UL)
            self.add_sound(audio_paths[5])
            self.play(Write(inflection_title))
            if real_inflection_points:
                inf_x = real_inflection_points[0]
                inf_y = cubic_function(inf_x)
                inflection_dot = Dot(axes.c2p(inf_x, inf_y), color=PURPLE)
                inflection_label = MathTex(f"({inf_x:.2f}, {inf_y:.2f})", font_size=24, color=PURPLE).next_to(inflection_dot, UP)
                slope_at_inflection = diff_function(inf_x)
                tangent_line = axes.plot(lambda x: inf_y + slope_at_inflection * (x - inf_x), x_range=[inf_x - 2, inf_x + 2], color=PURPLE)
                if a > 0:
                    concave_up = axes.plot(cubic_function, x_range=[inf_x, x_range[1]], color=RED)
                    concave_down = axes.plot(cubic_function, x_range=[x_range[0], inf_x], color=BLUE)
                else:
                    concave_down = axes.plot(cubic_function, x_range=[inf_x, x_range[1]], color=BLUE)
                    concave_up = axes.plot(cubic_function, x_range=[x_range[0], inf_x], color=RED)
                self.play(Create(inflection_dot), Write(inflection_label), Transform(graph, VGroup(concave_up, concave_down)))
                self.play(Create(tangent_line))
                concave_up_label = Text("Concave Up", font_size=20, color=RED).next_to(axes.c2p(inf_x + 1, cubic_function(inf_x + 1)), UR)
                concave_down_label = Text("Concave Down", font_size=20, color=BLUE).next_to(axes.c2p(inf_x - 1, cubic_function(inf_x - 1)), UL)
                self.play(Write(concave_up_label), Write(concave_down_label))
                self.wait(2)
                self.play(FadeOut(inflection_dot), FadeOut(inflection_label), FadeOut(tangent_line), FadeOut(concave_up_label), FadeOut(concave_down_label))
            else:
                inflection_message = Text("No inflection points exist", font_size=24, color=PURPLE).next_to(axes.get_center(), UP)
                self.play(Write(inflection_message))
                self.wait(1)
                self.play(FadeOut(inflection_message))
            self.play(FadeOut(inflection_title))
            
            # Rate of Change
            rate_title = Text("Rate of Change", font_size=32).to_corner(UL)
            self.add_sound(audio_paths[6])
            self.play(Write(rate_title))
            standard_graph = axes.plot(cubic_function, color=GREEN)
            self.play(Transform(graph, standard_graph))
            x_tracker = ValueTracker(x_range[0])
            moving_dot = always_redraw(lambda: Dot(axes.c2p(x_tracker.get_value(), cubic_function(x_tracker.get_value())), color=YELLOW))
            tangent_line = always_redraw(lambda: axes.plot(lambda t: cubic_function(x_tracker.get_value()) + diff_function(x_tracker.get_value()) * (t - x_tracker.get_value()), x_range=[x_tracker.get_value() - 1, x_tracker.get_value() + 1], color=YELLOW))
            slope_label = always_redraw(lambda: MathTex(f"\\text{{Slope}} = {diff_function(x_tracker.get_value()):.2f}", font_size=24, color=YELLOW).next_to(moving_dot, UP))
            self.play(Create(moving_dot), Create(tangent_line), Write(slope_label))
            self.play(x_tracker.animate.set_value(x_range[1]), run_time=6, rate_func=linear)
            self.wait(1)
            self.play(FadeOut(moving_dot), FadeOut(tangent_line), FadeOut(slope_label), FadeOut(rate_title))
            
            # Area Under the Curve
            area_title = Text("Area Under the Curve", font_size=32).to_corner(UL)
            self.play(Write(area_title))
            if len(real_roots) >= 2:
                sorted_roots = sorted(real_roots)
                for i in range(len(sorted_roots) - 1):
                    left_bound, right_bound = sorted_roots[i], sorted_roots[i + 1]
                    area = axes.get_area(graph, x_range=[left_bound, right_bound], color=BLUE if i % 2 == 0 else RED, opacity=0.5)
                    self.play(Create(area))
                    area_label = MathTex(f"\\int_{{{left_bound:.2f}}}^{{{right_bound:.2f}}} f(x) \\, dx", font_size=24).next_to(axes.c2p((left_bound + right_bound) / 2, -1), DOWN)
                    self.play(Write(area_label))
                    area_value = np.trapz([cubic_function(x) for x in np.linspace(left_bound, right_bound, 100)], np.linspace(left_bound, right_bound, 100))
                    area_value_label = MathTex(f"\\text{{Area}} = {abs(area_value):.2f}", font_size=24, color=BLUE if i % 2 == 0 else RED).next_to(area_label, DOWN)
                    self.play(Write(area_value_label))
                    self.wait(1)
                    self.play(FadeOut(area), FadeOut(area_label), FadeOut(area_value_label))
            else:
                area_message = Text("Not enough real roots to show a bounded area", font_size=24).next_to(axes.get_center(), DOWN)
                self.play(Write(area_message))
                self.wait(1)
                self.play(FadeOut(area_message))
            self.play(FadeOut(area_title))
            
            # Summary
            summary_title = Text("Summary", font_size=36).to_edge(UP)
            summary_points = [
                Text(f"• End behavior: {behavior}", font_size=24),
                Text(f"• Roots: {roots_str}", font_size=24),
                Text(f"• Turning points: {critical_str}", font_size=24),
                Text(f"• Inflection: {inflection_str}", font_size=24)
            ]
            summary_vgroup = VGroup(*summary_points).arrange(DOWN).next_to(summary_title, DOWN, buff=0.5)
            self.play(FadeOut(graph_group))
            self.play(FadeIn(summary_title))
            for point in summary_points:
                self.play(Write(point), run_time=0.5)
            self.wait(2)
            self.play(FadeOut(summary_title), FadeOut(summary_vgroup))
            self.wait(1)
    
    scene = CubicFunctionExplanation()
    scene.render()

if __name__ == "__main__":
    equation = "x^3-x^2+1"
    generate_video_cubic(equation)
    print("done")
