import streamlit as st
import requests

# Page Config
st.set_page_config(page_title="CuraBot - AI Medical Research", page_icon="⚕️")

st.title("⚕️ CuraBot")
st.subheader("AI-Powered Symptom Research Assistant")

# User Input
user_input = st.text_area("Describe your symptoms:", placeholder="e.g., persistent cough and mild fever")

if st.button("Analyze Symptoms"):
    if user_input:
        with st.spinner("🔍 Agents are searching PubMed..."):
            try:
                # This connects to your FastAPI backend
                response = requests.get(f"http://127.0.0.1:8000/analyze", params={"symptoms": user_input})
                
                if response.status_code == 200:
                    result = response.json()
                    st.success("Analysis Complete")
                    st.markdown("### 📋 Research Findings")
                    st.write(result.get("analysis", "No analysis found."))
                    st.warning("**Disclaimer:** For research only. Not a medical diagnosis.")
                else:
                    st.error(f"Backend Error: {response.status_code}")
            except Exception as e:
                st.error(f"Could not connect to backend. Is it running? Error: {e}")
    else:
        st.info("Please enter some symptoms.")