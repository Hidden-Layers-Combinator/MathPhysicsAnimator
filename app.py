# app.py
import streamlit as st
import os
from sympy import symbols, solve, sympify
from generate_video import generate_video

st.title("Quadratic Equation Graph Generator")

equation = st.text_input("Enter a quadratic equation (e.g., x^2 + 2x + 1):")
if st.button("Generate Video"):
    if equation:
      placeholder = st.empty()
      with placeholder:
        st.write("Generating video, please wait...")
      try:
        video_path,audio_path1,audio_path2,audio_path3,audio_path4= generate_video(equation)

            # Display video
        st.video(video_path)
            # Offer download
        with open(video_path, "rb") as f:
                st.download_button("Download Video", f, file_name="graph.mp4")
        with placeholder:
              st.write("Video generated successfully!")
            # Clean up
        os.remove(video_path)
        os.remove(audio_path1)
        os.remove(audio_path2)
        os.remove(audio_path3)
        os.remove(audio_path4)
      except ValueError as e:
        st.error(str(e))
    else:
      st.write("Please enter an equation.")

