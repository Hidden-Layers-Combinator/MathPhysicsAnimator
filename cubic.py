from manim import *
import numpy as np

class CubicFunctionExplanation(Scene):
    def construct(self):
        # Configuration
        axes_config = {
            "x_range": [-5, 5, 1],
            "y_range": [-8, 8, 2],
            "axis_config": {"color": BLUE},
            "x_axis_config": {"numbers_to_include": np.arange(-4, 5, 2)},
            "y_axis_config": {"numbers_to_include": np.arange(-8, 6, 2)},
        }
        
        # Create coordinate system
        axes = Axes(**axes_config)
        axes_labels = axes.get_axis_labels(x_label="x", y_label="f(x)").scale(0.5)
        
        # Initialize parameters for the cubic function
        a, b, c, d = 1, 0, -3, 0
        # Title
        title = Text("Understanding Cubic Equations", font_size=48)
        subtitle = Text(f"f(x) = {a}x³ + {b}x² + {c}x + {d}", font_size=36)
        subtitle.next_to(title, DOWN)
        
        self.play(AddTextLetterByLetter(title), run_time=1.5)
        self.play(AddTextLetterByLetter(subtitle), run_time=1)
        self.wait(0.5)

        self.play(
            title.animate.scale(1.2),
            subtitle.animate.scale(1.2),
            run_time=0.5
        )
        self.play(
            title.animate.scale(0.8333).to_edge(UP),  # 0.8333 = 1/1.2 to return to original size
            FadeOut(subtitle),
            FadeOut(title),
            run_time=0.7
        )
        
        # Display coordinate system
        self.play(Create(axes), Write(axes_labels))
        self.wait(1)
        
        # Initial cubic function
        def cubic_function(x, a=1, b=0, c=0, d=0):
            return a * x**3 + b * x**2 + c * x + d
        

        


        # Create function graph
        function_label = MathTex(f"f(x) = {a}x^3 {'' if b == 0 else '+ ' + str(b) + 'x^2'} {'' if c == 0 else ('- ' if c < 0 else '+ ') + str(abs(c)) + 'x'} {'' if d == 0 else ('- ' if d < 0 else '+ ') + str(abs(d))}")
        function_label.to_edge(UP)
        graph = axes.plot(lambda x: cubic_function(x, a, b, c, d), color=GREEN)

        
        # Display graph
        self.play(Create(graph))
        self.play(Write(function_label))
        self.wait(1)
        
        # 1. Shape and General Behavior Based on Leading Coefficient
        # shape_title = Text("Shape and General Behavior", font_size=32)
        # shape_title.to_corner(UL)
        # self.play(Write(shape_title))
        
        # Animate coefficient changing
        # explanation = Text("The leading coefficient determines end behavior", font_size=24)
        # explanation.next_to(shape_title, DOWN, aligned_edge=LEFT)
        # self.play(Write(explanation))
        
        # Animate changing the leading coefficient
        # for new_a in [2, -1, -2]:
        #     new_graph = axes.plot(lambda x: cubic_function(x, new_a, b, c, d), color=GREEN)
        #     new_function_label = MathTex(f"f(x) = {new_a}x^3 {'' if b == 0 else '+ ' + str(b) + 'x^2'} {'' if c == 0 else ('- ' if c < 0 else '+ ') + str(abs(c)) + 'x'} {'' if d == 0 else ('- ' if d < 0 else '+ ') + str(abs(d))}")
        #     new_function_label.next_to(title, DOWN)
            
        #     self.play(
        #         Transform(graph, new_graph),
        #         Transform(function_label, new_function_label)
        #     )
            
            # Create arrows showing behavior at infinities
        if a > 0:
            right_arrow = Arrow(axes.c2p(4, 7), axes.c2p(5, 8), color=YELLOW)
            left_arrow = Arrow(axes.c2p(-4, -7), axes.c2p(-5, -8), color=YELLOW)
            right_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
            left_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
        else:
            right_arrow = Arrow(axes.c2p(4, -7), axes.c2p(5, -8), color=YELLOW)
            left_arrow = Arrow(axes.c2p(-4, 7), axes.c2p(-5, 8), color=YELLOW)
            right_text = Text("→ -∞", font_size=20, color=YELLOW).next_to(right_arrow, RIGHT)
            left_text = Text("→ +∞", font_size=20, color=YELLOW).next_to(left_arrow, LEFT)
                
        self.play(
            Create(right_arrow),
            Create(left_arrow),
            Write(right_text),
            Write(left_text)
        )
        self.wait(1)
        self.play(
            FadeOut(right_arrow),
            FadeOut(left_arrow),
            FadeOut(right_text),
            FadeOut(left_text)
        )
        
        # Reset to standard form for next explanations
        a, b, c, d = 1, 0, -3, 0
        standard_graph = axes.plot(lambda x: cubic_function(x, a, b, c, d), color=GREEN)
        standard_function_label = MathTex(f"f(x) = {a}x^3 {'' if b == 0 else '+ ' + str(b) + 'x^2'} {'' if c == 0 else ('- ' if c < 0 else '+ ') + str(abs(c)) + 'x'} {'' if d == 0 else ('- ' if d < 0 else '+ ') + str(abs(d))}")
        standard_function_label.next_to(title, DOWN)
        
        self.play(
            Transform(graph, standard_graph),
            Transform(function_label, standard_function_label)
        )
        
        # 2. Number of Real Roots
        roots_title = Text("Real Roots", font_size=32)
        roots_title.to_corner(UL)
        self.play(Write(roots_title))
        
        # Find roots of the current function
        roots = [np.sqrt(3), 0, -np.sqrt(3)]  # Roots of x³-3x
        root_dots = []
        root_labels = []
        
        # Animate finding the roots
        for i, root in enumerate(roots):
            dot = Dot(axes.c2p(root, 0), color=RED)
            label = MathTex(f"x_{i+1} = {root:.2f}", color=RED, font_size=24)
            label.next_to(dot, DOWN)
            
            root_dots.append(dot)
            root_labels.append(label)
            
            self.play(
                Create(dot),
                Write(label)
            )
        
        self.wait(1)
        
        # 3. Turning Points (Local Maxima and Minima)
        self.play(
            *[FadeOut(label) for label in root_labels],
            *[FadeOut(dot) for dot in root_dots],
            FadeOut(roots_title)
        )
        
        turning_points_title = Text("Turning Points", font_size=32)
        turning_points_title.to_corner(UL)
        self.play(Write(turning_points_title))
        
        # Derivative function for slope
        def derivative(x):
            return 3*a*x**2 + 2*b*x + c
        
        # Find turning points (where derivative = 0)
        turning_points = [-1, 1]  # For x³-3x, derivative = 3x²-3 = 0 at x = ±1
        
        # Animate tangent lines at turning points
        for tp in turning_points:
            # Value at turning point
            y_tp = cubic_function(tp, a, b, c, d)
            
            # Create dot at turning point
            tp_dot = Dot(axes.c2p(tp, y_tp), color=YELLOW)
            
            # Label with coordinates
            tp_label = MathTex(f"({tp}, {y_tp:.2f})", color=YELLOW, font_size=24)
            tp_label.next_to(tp_dot, UP)
            
            # Create tangent line (slope = 0 at turning points)
            tangent_line = axes.plot(lambda x: y_tp, x_range=[tp-1, tp+1], color=YELLOW)
            
            self.play(
                Create(tp_dot),
                Write(tp_label)
            )
            
            # Moving tangent line animation before reaching turning point
            for x in np.linspace(tp-0.5, tp, 5):
                y = cubic_function(x, a, b, c, d)
                slope = derivative(x)
                
                # Line equation: y - y1 = m(x - x1)
                def tangent_at_x(t):
                    return y + slope * (t - x)
                
                temp_tangent = axes.plot(tangent_at_x, x_range=[x-1, x+1], color=BLUE)
                slope_label = MathTex(f"\\text{{Slope}} = {slope:.2f}", font_size=20, color=BLUE)
                slope_label.next_to(axes.c2p(x, y), UR)
                
                self.play(
                    Create(temp_tangent),
                    Write(slope_label),
                    run_time=0.3
                )
                self.play(
                    FadeOut(temp_tangent),
                    FadeOut(slope_label),
                    run_time=0.3
                )
            
            self.play(Create(tangent_line))
            self.wait(0.5)
            
            # Moving tangent line animation after turning point
            for x in np.linspace(tp, tp+0.5, 5):
                y = cubic_function(x, a, b, c, d)
                slope = derivative(x)
                
                # Line equation: y - y1 = m(x - x1)
                def tangent_at_x(t):
                    return y + slope * (t - x)
                
                temp_tangent = axes.plot(tangent_at_x, x_range=[x-1, x+1], color=BLUE)
                slope_label = MathTex(f"\\text{{Slope}} = {slope:.2f}", font_size=20, color=BLUE)
                slope_label.next_to(axes.c2p(x, y), UR)
                
                self.play(
                    Create(temp_tangent),
                    Write(slope_label),
                    run_time=0.3
                )
                self.play(
                    FadeOut(temp_tangent),
                    FadeOut(slope_label),
                    run_time=0.3
                )
            
            if tp == -1:
                max_min_label = Text("Local Maximum", font_size=20, color=YELLOW)
            else:
                max_min_label = Text("Local Minimum", font_size=20, color=YELLOW)
                
            max_min_label.next_to(tp_label, UP)
            self.play(Write(max_min_label))
            
            self.wait(1)
            self.play(
                FadeOut(tangent_line),
                FadeOut(tp_dot),
                FadeOut(tp_label),
                FadeOut(max_min_label)
            )
        
        self.play(FadeOut(turning_points_title))
        
        # 4. Inflection Point
        inflection_title = Text("Inflection Point", font_size=32)
        inflection_title.to_corner(UL)
        self.play(Write(inflection_title))
        
        # Inflection point (where second derivative = 0)
        inflection_x = 0  # For x³-3x, second derivative = 6x = 0 at x = 0
        inflection_y = cubic_function(inflection_x, a, b, c, d)
        
        # Create dot at inflection point
        inflection_dot = Dot(axes.c2p(inflection_x, inflection_y), color=PURPLE)
        inflection_label = MathTex(f"\\text{{Inflection at }} ({inflection_x}, {inflection_y})", font_size=24, color=PURPLE)
        inflection_label.next_to(inflection_dot, UP)
        
        # Show different colored segments for different concavity
        concave_up = axes.plot(lambda x: cubic_function(x, a, b, c, d), x_range=[0, 5], color=RED)
        concave_down = axes.plot(lambda x: cubic_function(x, a, b, c, d), x_range=[-5, 0], color=BLUE)
        
        self.play(
            Create(inflection_dot),
            Write(inflection_label),
            Transform(graph, VGroup(concave_up, concave_down))
        )
        
        # Tangent line at inflection point
        slope_at_inflection = derivative(inflection_x)
        
        def tangent_at_inflection(x):
            return inflection_y + slope_at_inflection * (x - inflection_x)
            
        tangent_line = axes.plot(tangent_at_inflection, x_range=[-2, 2], color=PURPLE)
        
        self.play(Create(tangent_line))
        
        # Label concavity
        concave_up_label = Text("Concave Up", font_size=20, color=RED)
        concave_up_label.next_to(axes.c2p(2, cubic_function(2, a, b, c, d)), UR)
        
        concave_down_label = Text("Concave Down", font_size=20, color=BLUE)
        concave_down_label.next_to(axes.c2p(-2, cubic_function(-2, a, b, c, d)), UL)
        
        self.play(
            Write(concave_up_label),
            Write(concave_down_label)
        )
        
        self.wait(2)
        
        self.play(
            FadeOut(inflection_dot),
            FadeOut(inflection_label),
            FadeOut(tangent_line),
            FadeOut(concave_up_label),
            FadeOut(concave_down_label),
            FadeOut(inflection_title)
        )
        
        # 5. Rate of Change Animation
        rate_title = Text("Rate of Change", font_size=32)
        rate_title.to_corner(UL)
        self.play(Write(rate_title))
        
        # Reset the graph
        standard_graph = axes.plot(lambda x: cubic_function(x, a, b, c, d), color=GREEN)
        self.play(Transform(graph, standard_graph))
        
        # Create a moving point and tangent line
        x_tracker = ValueTracker(-3)
        
        def get_x():
            return x_tracker.get_value()
            
        def get_y():
            return cubic_function(get_x(), a, b, c, d)
            
        def get_slope():
            return derivative(get_x())
        
        # Moving dot
        moving_dot = always_redraw(
            lambda: Dot(axes.c2p(get_x(), get_y()), color=YELLOW)
        )
        
        # Function for tangent line
        def get_tangent_line():
            x = get_x()
            y = get_y()
            slope = get_slope()
            
            def tangent_func(t):
                return y + slope * (t - x)
                
            return axes.plot(tangent_func, x_range=[x-1, x+1], color=YELLOW)
            
        tangent_line = always_redraw(get_tangent_line)
        
        # Slope label
        slope_label = always_redraw(
            lambda: MathTex(
                f"\\text{{Slope}} = {get_slope():.2f}",
                font_size=24,
                color=YELLOW
            ).next_to(moving_dot, UP)
        )
        
        self.play(
            Create(moving_dot),
            Create(tangent_line),
            Write(slope_label)
        )
        
        # Animate the point moving along the curve
        self.play(
            x_tracker.animate.set_value(3),
            run_time=6,
            rate_func=linear
        )
        
        self.wait(1)
        
        self.play(
            FadeOut(moving_dot),
            FadeOut(tangent_line),
            FadeOut(slope_label),
            FadeOut(rate_title)
        )
        
        # 6. Area Under the Curve
        area_title = Text("Area Under the Curve", font_size=32)
        area_title.to_corner(UL)
        self.play(Write(area_title))
        
        # Reset function to show area between roots
        # Shade area between roots
        area = axes.get_area(graph, x_range=[-np.sqrt(3), 0], color=BLUE, opacity=0.5)
        area2 = axes.get_area(graph, x_range=[0, np.sqrt(3)], color=RED, opacity=0.5)
        
        # Animate shading
        self.play(Create(area))
        area_label = MathTex("\\int_{-\\sqrt{3}}^{0} (x^3-3x) \\, dx", font_size=24)
        area_label.next_to(axes.c2p(-1, -1), DOWN)
        self.play(Write(area_label))
        
        self.play(Create(area2))
        area_label2 = MathTex("\\int_{0}^{\\sqrt{3}} (x^3-3x) \\, dx", font_size=24)
        area_label2.next_to(axes.c2p(1, -1), DOWN)
        self.play(Write(area_label2))
        
        # Calculate and show the area
        # For x³-3x, integral is x⁴/4 - 3x²/2
        def integral(x):
            return x**4/4 - 3*x**2/2
            
        area_value = integral(0) - integral(-np.sqrt(3))
        area_value2 = integral(np.sqrt(3)) - integral(0)
        
        area_value_label = MathTex(f"\\text{{Area}} = {abs(area_value):.2f}", font_size=24, color=BLUE)
        area_value_label.next_to(area_label, DOWN)
        
        area_value_label2 = MathTex(f"\\text{{Area}} = {abs(area_value2):.2f}", font_size=24, color=RED)
        area_value_label2.next_to(area_label2, DOWN)
        
        self.play(
            Write(area_value_label),
            Write(area_value_label2)
        )
        
        self.wait(2)
        
        # Final summary
        self.play(
            FadeOut(area),
            FadeOut(area2),
            FadeOut(area_label),
            FadeOut(area_label2),
            FadeOut(area_value_label),
            FadeOut(area_value_label2),
            FadeOut(area_title)
        )
        
        summary_title = Text("Summary of Cubic Functions", font_size=36)
        summary_title.to_edge(UP)
        
        summary_points = VGroup(
            Text("• Shape depends on leading coefficient", font_size=24),
            Text("• Up to 3 real roots", font_size=24),
            Text("• 2 turning points (max and min)", font_size=24),
            Text("• 1 inflection point", font_size=24),
            Text("• Cubic growth rate", font_size=24)
        ).arrange(DOWN, aligned_edge=LEFT)
        summary_points.next_to(summary_title, DOWN, buff=0.5)
        
        self.play(
            FadeOut(title),
            FadeOut(function_label),
            FadeIn(summary_title),
        )
        
        for point in summary_points:
            self.play(Write(point), run_time=0.5)
            
        self.wait(2)
        
        # Final fade out
        self.play(
            FadeOut(summary_title),
            FadeOut(summary_points),
            FadeOut(axes),
            FadeOut(axes_labels),
            FadeOut(graph)
        )
        
        self.wait(1)
scene = CubicFunctionExplanation()
scene.render()