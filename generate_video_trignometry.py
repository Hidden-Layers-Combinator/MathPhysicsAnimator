from manim import *
import numpy as np
from sympy import *
import re

class TrigFunctionVisualizer(Scene):
    def construct(self):
        # This function can be updated with the user's equation
        def get_user_equation():
            # Replace this with actual user input
            # Example: "2*sin(x) + cos(2*x)"
            return "2*sin(x) + cos(2*x)"
        
        # Parse and process the user's equation
        equation_str = get_user_equation()
        
        # Process the equation with sympy for analysis
        x = Symbol('x')
        expr = sympify(equation_str)
        
        # Function to convert sympy expression to a manim-compatible function
        def expr_to_func(expression):
            return lambda x_val: float(expression.subs(x, x_val))
        
        func = expr_to_func(expr)
        
        # Calculate important properties of the function
        def find_period():
            # Simplified approach for period finding
            if "sin" in equation_str or "cos" in equation_str:
                coefficients = []
                for term in ["sin", "cos"]:
                    # Pattern for forms like sin(2x) or sin(2*x)
                    pattern = fr"{term}\((\d*\.?\d*)[\*]?x\)"
                    matches = re.findall(pattern, equation_str)
                    if matches:
                        for match in matches:
                            if match == '':
                                coefficients.append(1.0)
                            else:
                                coefficients.append(float(match))
                    # Pattern for forms like sin(x*2) or sin(x*2.5)
                    pattern = fr"{term}\(x[\*]?(\d*\.?\d*)\)"
                    matches = re.findall(pattern, equation_str)
                    if matches:
                        for match in matches:
                            if match == '':
                                coefficients.append(1.0)
                            else:
                                coefficients.append(float(match))
                
                if coefficients:
                    if len(coefficients) > 1:
                        from math import gcd
                        from functools import reduce
                        from fractions import Fraction
                        
                        # Convert coefficients to fractions to handle floats
                        def lcm(a, b):
                            # Convert to integers via fractions
                            fa, fb = Fraction(a).limit_denominator(), Fraction(b).limit_denominator()
                            num = abs(fa.numerator * fb.numerator)
                            denom = gcd(fa.denominator, fb.denominator)
                            lcm_num = num // gcd(num, denom)
                            return lcm_num / (fa.denominator * fb.denominator // denom)
                        
                        def lcm_multiple(numbers):
                            return reduce(lcm, numbers)
                        
                        coefficient = lcm_multiple(coefficients)
                    else:
                        coefficient = coefficients[0]
                    
                    return 2 * PI / coefficient
                return 2 * PI
            elif "tan" in equation_str:
                return PI
            return 2 * PI  # Default
        
        period = find_period()
        
        # For analysis, evaluate over a range
        x_vals = np.linspace(0, 2 * period, 1000)
        y_vals = [func(x_val) for x_val in x_vals]
        
        # Find amplitude/range
        min_y = min(y_vals)
        max_y = max(y_vals)
        amplitude = (max_y - min_y) / 2

        # Find zeros (approximate)
        zeros = []
        for i in range(1, len(x_vals)):
            if (y_vals[i-1] <= 0 and y_vals[i] >= 0) or (y_vals[i-1] >= 0 and y_vals[i] <= 0):
                t = -y_vals[i-1] / (y_vals[i] - y_vals[i-1])
                zero_x = x_vals[i-1] + t * (x_vals[i] - x_vals[i-1])
                zeros.append(zero_x)
        
        # Find local maxima and minima (approximate)
        critical_points = []
        for i in range(1, len(y_vals) - 1):
            if (y_vals[i-1] < y_vals[i] and y_vals[i] > y_vals[i+1]) or \
               (y_vals[i-1] > y_vals[i] and y_vals[i] < y_vals[i+1]):
                critical_points.append((x_vals[i], y_vals[i]))
        
        # Setup the visualization
        x_min, x_max = 0, 2 * period
        y_padding = max(1, (max_y - min_y) * 0.2)
        y_min, y_max = min_y - y_padding, max_y + y_padding
        
        axes = Axes(
            x_range=[x_min, x_max, PI/2],
            y_range=[y_min, y_max, max(1, (max_y - min_y) / 4)],
            axis_config={"color": BLUE},
            # x_axis_config={
            #     "numbers_to_include": np.arange(0, x_max + 0.1, PI/2),
            #     "numbers_with_elongated_ticks": np.arange(0, x_max + 0.1, PI)
            # }
        )
        
        # Custom x-axis labels
        x_labels = []
        for i in range(int(x_max / (PI/2)) + 1):
            value = i * PI/2
            if i == 0:
                label = MathTex("0")
            elif i % 4 == 0:
                label = MathTex(fr"{i // 2}\pi")
            elif i % 2 == 0:
                label = MathTex(fr"{i // 2}\pi")
            else:
                label = MathTex(fr"\frac{{{i}}}{{2}}\pi")
            label.scale(0.5).next_to(axes.c2p(value, 0), DOWN)
            x_labels.append(label)
        
        axes_labels = axes.get_axis_labels(
            x_label=MathTex("x"),
            y_label=MathTex("f(x)")
        )
        
        # Title and equation display
        title = Title(f"Trigonometric Function Analysis", include_underline=False)
        equation_display = MathTex(f"f(x) = {equation_str}")
        equation_display.next_to(title, DOWN)
        
        # Create the graph
        graph = axes.plot(func, x_range=[x_min, x_max], color=WHITE)
        graph_group=VGroup(graph,axes,axes_labels,x_labels).shift(DOWN*0.7)
        
        # Initial scene setup
        self.play(
            Write(title),
            Write(equation_display)
        )
        self.play(
            Create(axes),
            *[Write(label) for label in x_labels],
            Write(axes_labels)
        )
        self.wait(1)
        
        # Create the graph with a growing animation
        self.play(Create(graph), run_time=2)
        self.wait(1)
        self.play(
            FadeOut(title),
            FadeOut(equation_display)
        )
        # 1. Periodicity
        periodicity_title = Text("Periodicity", font_size=30).to_edge(UP)
        self.play(Write(periodicity_title))
        
        period_segment = axes.plot(func, x_range=[0, period], color=YELLOW)
        period_label = MathTex(fr"\text{{Period}} = {period/PI:.2f}\pi").next_to(periodicity_title, DOWN, aligned_edge=LEFT)
        
        self.play(Create(period_segment))
        self.play(Write(period_label))
        
        second_period = axes.plot(func, x_range=[period, 2*period], color=GREEN)
        self.play(Create(second_period))
        
        start_point = axes.c2p(0, func(0))
        mid_point = axes.c2p(period, func(period))
        end_point = axes.c2p(2*period, func(2*period))
        
        period_arrow1 = Arrow(start_point, mid_point, color=YELLOW, buff=0.1)
        period_arrow2 = Arrow(mid_point, end_point, color=GREEN, buff=0.1)
        
        self.play(Create(period_arrow1))
        self.play(Create(period_arrow2))
        
        self.wait(1)
        
        self.play(
            FadeOut(period_segment),
            FadeOut(period_arrow1),
            FadeOut(period_arrow2),
            FadeOut(second_period),
            FadeOut(periodicity_title),
            FadeOut(period_label)
        )
        
        # 2. Amplitude and Range
        amplitude_title = Text("Amplitude and Range", font_size=30).to_edge(UP)
        self.play(Write(amplitude_title))
        
        range_text = MathTex(fr"\text{{Range}} = [{min_y:.2f}, {max_y:.2f}]")
        amplitude_text = MathTex(fr"\text{{Amplitude}} \approx {amplitude:.2f}")
        
        range_text.next_to(amplitude_title, DOWN, aligned_edge=LEFT)
        amplitude_text.next_to(range_text, DOWN, aligned_edge=LEFT)
        
        max_x = x_vals[np.argmax(y_vals)]
        min_x = x_vals[np.argmin(y_vals)]
        
        max_point = Dot(axes.c2p(max_x, max_y), color=RED)
        min_point = Dot(axes.c2p(min_x, min_y), color=BLUE)
        
        range_arrow = DoubleArrow(
            axes.c2p(max_x, min_y),
            axes.c2p(max_x, max_y),
            color=YELLOW,
            buff=0.1
        )
        
        self.play(Write(range_text), Write(amplitude_text))
        self.play(Create(max_point), Create(min_point))
        self.play(Create(range_arrow))
        
        self.wait(1)
        
        self.play(
            FadeOut(amplitude_title),
            FadeOut(range_text),
            FadeOut(amplitude_text),
            FadeOut(max_point),
            FadeOut(min_point),
            FadeOut(range_arrow)
        )
        
        # 3. Zeros (X-Intercepts)
        zeros_title = Text("Zeros (X-Intercepts)", font_size=36).to_corner(UL)
        self.play(Write(zeros_title))
        
        zero_dots = []
        zero_labels = []
        
        visible_zeros = [z for z in zeros if x_min <= z <= x_max]
        max_zeros_to_show = 5
        if len(visible_zeros) > max_zeros_to_show:
            visible_zeros = visible_zeros[:max_zeros_to_show]
        
        for i, zero in enumerate(visible_zeros):
            dot = Dot(axes.c2p(zero, 0), color=RED)
            zero_dots.append(dot)
            
            pi_multiple = zero / PI
            if abs(pi_multiple - round(pi_multiple)) < 0.01:
                if round(pi_multiple) == 0:
                    label_text = "x = 0"
                elif round(pi_multiple) == 1:
                    label_text = r"x = \pi"
                else:
                    label_text = fr"x = {round(pi_multiple)}\pi"
            else:
                denominator = 2
                while denominator <= 6:
                    if abs(pi_multiple * denominator - round(pi_multiple * denominator)) < 0.05:
                        numerator = round(pi_multiple * denominator)
                        if numerator == 0:
                            label_text = "x = 0"
                        elif numerator == denominator:
                            label_text = r"x = \pi"
                        else:
                            label_text = fr"x = \frac{{{numerator}}}{{{denominator}}}\pi"
                        break
                    denominator += 1
                else:
                    label_text = f"x = {zero:.2f}"
            
            label = MathTex(label_text)
            label.scale(0.7).next_to(dot, DOWN)
            zero_labels.append(label)
            
            self.play(
                Create(dot),
                Write(label)
            )
        
        self.wait(1)
        
        self.play(
            FadeOut(zeros_title),
            *[FadeOut(dot) for dot in zero_dots],
            *[FadeOut(label) for label in zero_labels]
        )
        
        # 4. Maxima and Minima
        extrema_title = Text("Maxima and Minima", font_size=36).to_corner(UL)
        self.play(Write(extrema_title))
        
        max_cp_to_show = 6
        if len(critical_points) > max_cp_to_show:
            critical_points = critical_points[:max_cp_to_show]
        
        for i, (cp_x, cp_y) in enumerate(critical_points):
            cp_dot = Dot(axes.c2p(cp_x, cp_y), color=YELLOW)
            tangent = axes.plot(
                lambda x: cp_y,
                x_range=[cp_x - 0.2 * period, cp_x + 0.2 * period],
                color=GREEN
            )
            
            cp_type = "Maximum" if i % 2 == 0 else "Minimum"
            pi_multiple = cp_x / PI
            if abs(pi_multiple - round(pi_multiple)) < 0.01:
                if round(pi_multiple) == 0:
                    x_text = "0"
                elif round(pi_multiple) == 1:
                    x_text = r"\pi"
                else:
                    x_text = fr"{round(pi_multiple)}\pi"
            else:
                denominator = 2
                while denominator <= 6:
                    if abs(pi_multiple * denominator - round(pi_multiple * denominator)) < 0.05:
                        numerator = round(pi_multiple * denominator)
                        if numerator == 0:
                            x_text = "0"
                        elif numerator == denominator:
                            x_text = r"\pi"
                        else:
                            x_text = fr"\frac{{{numerator}}}{{{denominator}}}\pi"
                        break
                    denominator += 1
                else:
                    x_text = f"{cp_x:.2f}"
            
            cp_label = MathTex(fr"{cp_type} \ at \ x = {x_text}")
            value_label = MathTex(fr"f({x_text}) = {cp_y:.2f}")
            
            cp_label.scale(0.6).next_to(cp_dot, UP)
            value_label.scale(0.6).next_to(cp_label, UP)
            
            self.play(
                Create(cp_dot),
                Create(tangent)
            )
            self.play(
                Write(cp_label),
                Write(value_label)
            )
            
            self.wait(0.5)
            
            self.play(
                FadeOut(cp_dot),
                FadeOut(tangent),
                FadeOut(cp_label),
                FadeOut(value_label)
            )
        
        self.play(FadeOut(extrema_title))
        
        # 5. Symmetry
        symmetry_title = Text("Symmetry", font_size=36).to_corner(UL)
        self.play(Write(symmetry_title))
        
        is_even = True
        is_odd = True
        test_points = np.linspace(0, period/2, 20)
        
        for t in test_points:
            if abs(func(t) - func(-t)) > 0.1:
                is_even = False
            if abs(func(t) + func(-t)) > 0.1:
                is_odd = False
                
        if is_even:
            symmetry_text = "Function has even symmetry: f(-x) = f(x)"
            reflect_point = 0
        elif is_odd:
            symmetry_text = "Function has odd symmetry: f(-x) = -f(x)"
            reflect_point = 0
        else:
            symmetry_text = "Function has no odd/even symmetry about origin"
            half_period_symmetry = True
            for t in test_points:
                if abs(func(t) - func(period/2 - t)) > 0.1:
                    half_period_symmetry = False
                    break
                    
            if half_period_symmetry:
                symmetry_text += fr"\nBut has symmetry about x = {period/(2*PI):.2f}\pi"
                reflect_point = period/2
            else:
                reflect_point = None
        
        symmetry_label = Text(symmetry_text, font_size=24)
        symmetry_label.next_to(symmetry_title, DOWN, aligned_edge=LEFT)
        
        self.play(Write(symmetry_label))
        
        if reflect_point is not None:
            if reflect_point == 0:
                if is_even:
                    original_segment = axes.plot(func, x_range=[0, period/2], color=YELLOW)
                    reflected_segment = axes.plot(lambda x: func(-x), x_range=[-period/2, 0], color=GREEN)
                else:
                    original_segment = axes.plot(func, x_range=[0, period/2], color=YELLOW)
                    reflected_segment = axes.plot(lambda x: -func(-x), x_range=[-period/2, 0], color=GREEN)
            else:
                original_segment = axes.plot(func, x_range=[0, reflect_point], color=YELLOW)
                reflected_segment = axes.plot(lambda x: func(2*reflect_point - x), 
                                              x_range=[reflect_point, 2*reflect_point], color=GREEN)
            
            self.play(Create(original_segment))
            self.play(Create(reflected_segment))
            
            self.wait(1)
            
            self.play(
                FadeOut(original_segment),
                FadeOut(reflected_segment)
            )
        
        self.wait(1)
        
        self.play(
            FadeOut(symmetry_title),
            FadeOut(symmetry_label)
        )
        
        # Final summary
        summary_title = Text("Summary", font_size=36).to_edge(UP)
        summary_points = VGroup(
            Text(f"• Period: {period/PI:.2f}π", font_size=24),
            Text(f"• Range: [{min_y:.2f}, {max_y:.2f}]", font_size=24),
            Text(f"• Amplitude: {amplitude:.2f}", font_size=24),
            Text(f"• {len(zeros)} zeros (x-intercepts) in [0, {2*period/PI:.0f}π]", font_size=24),
            Text(f"• {len(critical_points)} critical points", font_size=24),
            Text(f"• {'Has' if is_even or is_odd else 'No'} {'even' if is_even else 'odd' if is_odd else ''} symmetry", font_size=24)
        ).arrange(DOWN, aligned_edge=LEFT)
        
        summary_points.next_to(summary_title, DOWN, buff=0.5)
        
        self.play(
            FadeOut(title),
            FadeOut(equation_display),
            FadeIn(summary_title)
        )
        
        for point in summary_points:
            self.play(Write(point), run_time=0.5)
            
        self.wait(2)
        
        self.play(
            FadeOut(summary_title),
            FadeOut(summary_points),
            FadeOut(axes),
            # *[FadeOut(label) for label in x_labels],
            FadeOut(axes_labels),
            FadeOut(graph)
        )
        
        self.wait(1)

scene = TrigFunctionVisualizer()
scene.render()
