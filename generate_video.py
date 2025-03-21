from gtts import gTTS
from sympy import symbols, solve, sympify
from manim import *
from pydub import AudioSegment
import time
def generate_video(equation):
  x=symbols('x')
  eq=sympify(equation)
  a,b,c=eq.as_poly(x).all_coeffs()
  if a == 0:
        raise ValueError("Not a quadratic equation (a cannot be 0).")
  vertex_x = -b / (2 * a)
  vertex_y = a * vertex_x**2 + b * vertex_x + c
  discriminant = b**2 - 4 * a * c
  roots = solve(eq,x) if discriminant >= 0 else []
  #script
  wordEquation=''
  if a==1:
    wordEquation+="x square "
  elif a==-1:
    wordEquation+="minus x square "
  elif a<0:
    wordEquation+="minus "+str(abs(a))+" x square "
  else:
    wordEquation+=str(a)+" x square "
  if b==1:
    wordEquation+="plus x"
  elif b==-1:
    wordEquation+="minus x"
  elif b>0:
    wordEquation+="plus "+str(b)+" x"
  elif b<0:
    wordEquation+="minus "+str(abs(b))+" x"
  if c>0:
    wordEquation+=" plus "+str(c)
  elif c<0:
    wordEquation+=" minus "+str(abs(c))
  Intro=f'Welcome to Animath. I am your instructor jojo. You have given the equation {wordEquation}. Lets first note the coefficients '
  roots_str = (
        f" It has two real roots at x equals {roots[0]:.1f} and {roots[1]:.1f}." if len(roots) == 2 else
        f" It has one real root at x equals {roots[0]:.1f}." if len(roots) == 1 else
        " It has no real roots."
    )
  given=f'Here, a is equals {a}, b is equals {b} and c is equals {c}. '
  vertex=f'The vertex of the graph is at {vertex_x:.1f},{vertex_y:.1f}. '
  d=f'The discriminant is {discriminant}. '
  draw='Lets Draw the graph.'
  script=[Intro,given,vertex+d+roots_str]

  unique_id = time.time()
  audio_path1 = f"/tmp/audio1_{unique_id}.mp3"
  audio_path2 = f"/tmp/audio2_{unique_id}.mp3"
  audio_path3 = f"/tmp/audio3_{unique_id}.mp3"
  audio_path4 = f"/tmp/audio4_{unique_id}.mp3"
  #first audio
  first = gTTS(text=script[0], lang='en')
  first.save(audio_path1)
  #second audio
  second = gTTS(text=script[1], lang='en')
  second.save(audio_path2)
  #third audio
  third=gTTS(text=script[2],lang='en')
  third.save(audio_path3)
  #fourth audio
  fourth=gTTS(text=draw,lang='en')
  fourth.save(audio_path4)
  audio_1 = AudioSegment.from_mp3(audio_path1)
  audio_2=AudioSegment.from_mp3(audio_path2)
  audio_3=AudioSegment.from_mp3(audio_path3)
  audio_4=AudioSegment.from_mp3(audio_path4)

  # Combine audio with silence to match animation timing
  combined_audio = audio_1 + AudioSegment.silent(duration=1000)  # 1s pause after intro
  combined_audio += audio_2 
  combined_audio += audio_4+AudioSegment.silent(duration=3500)
  combined_audio += audio_3
  
  combined_audio_path = f"/tmp/combined_{unique_id}.mp3"
  combined_audio.export(combined_audio_path, format="mp3")
  class QuadraticGraph(MovingCameraScene):
    def construct(self):
      audio_length_1=len(audio_1)/1000
      audio_length_2=len(audio_2)/1000
      audio_length_3=len(audio_3)/1000
      audio_length_4=len(audio_4)/1000
      #Intro
      intro_text=Text(f'Welcome to Animath. I am your instructor jojo. You have given the equation {wordEquation}. Lets first note the coefficients').scale(0.3)
      #Given of Equation
      a_text = Text(f"a = {a}", color=WHITE).to_corner(UP + LEFT, buff=0.5)
      b_text = Text(f"b = {b}", color=WHITE).next_to(a_text, DOWN, aligned_edge=LEFT, buff=0.6)
      c_text = Text(f"c = {c}", color=WHITE).next_to(b_text, DOWN, aligned_edge=LEFT,buff=0.6)

      #Axes & Graph of equation
      axes = Axes(x_range=[-5, 5, 1], y_range=[-5, 5, 1], x_length=6, y_length=6)
      x_label = Text("x", font_size=24).next_to(axes.x_axis, RIGHT)
      y_label = Text("f(x)", font_size=24).next_to(axes.y_axis, UP)

      graph = axes.plot(lambda x: a * x**2 + b * x + c, color=BLUE)

      #grouping
      graph_group = VGroup(axes, x_label, y_label, graph, a_text, b_text, c_text)

      #vertex
      vertex = Dot(axes.coords_to_point(vertex_x, vertex_y), color=GREEN)
      vertex_label = Text(f"Vertex ({vertex_x},{vertex_y})").scale(0.3)
      vertex_label.next_to(vertex, DOWN)
      d_text=Text(f"Discriminant = {discriminant}").scale(0.6).to_corner(DOWN+RIGHT,buff=0.6)
      #roots
      if(len(roots)==1):
        root1=Dot(axes.coords_to_point(roots[0],0),color=YELLOW)
        root1_text=Text(f"({roots[0]},{0})").next_to(root1,UP).scale(0.4)
      elif(len(roots)==2):
        root1=Dot(axes.coords_to_point(roots[0],0),color=YELLOW)
        root1_text=Text(f"({roots[0]},{0})").next_to(root1,UP).scale(0.4)
        root2=Dot(axes.coords_to_point(roots[1],0),color=YELLOW)
        root2_text=Text(f"({roots[1]},{0})").next_to(root2,UP).scale(0.4)
      else:
        r_text=Text(f"Roots = {roots}").scale(0.6).next_to(UP,buff=0.6)
    
      #discriminant

      #Animation
      self.add_sound(combined_audio_path)
      self.play(AddTextLetterByLetter(intro_text),run_time=audio_length_1)
      self.play(FadeOut(intro_text))
      # self.add_sound(audio_path2)
      self.play(Write(a_text), Write(b_text), Write(c_text), run_time=audio_length_2)
      # self.add_sound(audio_path4)
      self.play(Create(axes), Write(x_label),Write(y_label),run_time=2)
      self.play(Create(graph),run_time=3)
      # self.add_sound(audio_path3)
      self.play(
            self.camera.frame.animate.set_width(2).move_to(vertex.get_center()),
            run_time=1
      )
      self.play(Create(vertex),run_time=1)
      self.play(Write(vertex_label),run_time=1)
      self.wait(2)
      self.play(
            self.camera.frame.animate.set_width(graph_group.width*1.5).move_to(axes.get_center()),run_time=0.5
        )
      self.play(Write(d_text),run_time=2)
      if (len(roots)==1):
        self.play(
            self.camera.frame.animate.set_width(1).move_to(root1.get_center()),
            run_time=0.5
        )
        self.play(

              Create(root1),
              Write(root1_text),
              run_time=1
          )
        self.wait()
        self.play(
              self.camera.frame.animate.set_width(graph_group.width * 1.5).move_to(axes.get_center()),
              run_time=0.5
          )
      elif(len(roots)==2):
        self.play(
            self.camera.frame.animate.set_width(1).move_to(root1.get_center()),
            run_time=0.5
        )
        self.play(
              Create(root1),
              Write(root1_text),
              run_time=1
          )
        self.play(
              self.camera.frame.animate.set_width(graph_group.width * 1.5).move_to(axes.get_center()),
              run_time=0.5
          )
        self.play(
            self.camera.frame.animate.set_width(1).move_to(root2.get_center()),
            run_time=0.5
        )
        self.play(
              Create(root2),
              Write(root2_text),
              run_time=1
          )
        self.play(
              self.camera.frame.animate.set_width(graph_group.width * 1.5).move_to(axes.get_center()),
              run_time=0.5
          )
      else:
        self.play(Write(r_text))
      self.wait(3)
  video_path = "/tmp/QuadraticGraph.mp4"
  config["output_file"] = video_path  # Tell Manim where to save the video
  scene = QuadraticGraph()
  scene.render(preview=False)
  return video_path, audio_path1, audio_path2, audio_path3, audio_path4




