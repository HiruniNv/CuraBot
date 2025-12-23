import streamlit as st
import requests

# Page Config
st.set_page_config(page_title="CuraBot - AI Medical Research", page_icon="⚕️", layout="centered")

# Custom CSS for a clean look
st.markdown("""
    <style>
    .main { background-color: #f5f7f9; }
    .stButton>button { width: 100%; border-radius: 5px; height: 3em; background-color: #007bff; color: white; }
    </style>
    """, unsafe_allow_html=True)

st.title("⚕️ CuraBot")
st.subheader("AI-Powered Symptom Research Assistant")
st.write("Describe how you're feeling, and our AI agents will search PubMed for related medical literature.")

# User Input
user_input = st.text_area("How are you feeling today?", placeholder="e.g., I've had a sharp pain in my lower back and a mild fever for two days.")

if st.button("Analyze Symptoms"):
    if user_input:
        with st.spinner("🔍 Agents are analyzing and searching PubMed..."):
            try:
                # Call your FastAPI Backend
                response = requests.get(f"http://127.0.0.1:8000/analyze", params={"symptoms": user_input})
                
                if response.status_code == 200:
                    result = response.json()
                    
                    st.success("Analysis Complete")
                    st.markdown("### 📋 Research Findings")
                    st.write(result.get("analysis", "No analysis found."))
                    
                    st.warning("**Disclaimer:** This information is for research purposes only and is not a medical diagnosis. Please consult a healthcare professional.")
                else:
                    st.error(f"Error: {response.status_code}")
            except Exception as e:
                st.error(f"Could not connect to backend: {e}")
    else:
        st.info("Please enter some symptoms first.")

st.divider()
st.caption("Powered by CrewAI, Groq (Llama 3.3), and PubMed Data.")