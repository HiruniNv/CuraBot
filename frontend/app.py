import streamlit as st
import requests
import pandas as pd
from fpdf import FPDF
from datetime import datetime
import time
import re

# --- 1. PAGE CONFIGURATION & SESSION STATE ---
# Sets up the browser tab title, favicon, and expands the layout to use full screen width.
st.set_page_config(
    page_title="CuraBot | AI Medical Research", 
    page_icon="⚕️", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize Session States: Manages user status, premium access, and app loading states.
states = {
    'is_logged_in': False,
    'is_premium': False,
    'history': [],
    'loading': False,
    'report_content': None
}
for key, val in states.items():
    if key not in st.session_state:
        st.session_state[key] = val

# --- 2. HELPER FUNCTIONS ---
# Validates email format using Regex to ensure a standard "user@domain.com" structure.
def is_valid_email(email):
    return re.match(r"[^@]+@[^@]+\.[^@]+", email)

# --- 3. CUSTOM CSS STYLING ---
# Injects custom CSS to define the dark medical theme and premium component styles.
cancel_color = "#dc3545" if st.session_state.loading else "#2c313c"

st.markdown(f"""
    <style>
    /* Dark background for the main application area */
    .main {{ background-color: #0e1117; }}
    
    /* Branding and Typography */
    .medical-title {{ font-size: 2.8rem; font-weight: 800; color: #ffffff; margin-bottom: 0px; }}
    .medical-subtitle {{ color: #4b9fff; font-size: 1.1rem; font-weight: 500; margin-bottom: 25px; }}
    
    /* Analyze Button: Primary Green Action */
    div.stButton > button:first-child[kind="primary"] {{
        background-color: #28a745 !important;
        color: white !important;
        border: none !important;
        height: 3.5em !important;
        width: 100% !important;
        font-weight: bold !important;
    }}
    
    /* Premium Feature list items inside the Go Premium popover */
    .premium-feature {{
        padding: 10px;
        border-radius: 8px;
        background: #1e2129;
        margin-bottom: 10px;
        border-left: 3px solid #d4af37;
    }}

    /* Container for AI-generated clinical findings */
    .result-card {{
        background-color: #1e2129;
        padding: 25px;
        border-radius: 15px;
        border-left: 5px solid #4b9fff;
        margin-top: 20px;
        line-height: 1.6;
    }}

    /* Red-tinted box for medical disclaimers and warnings */
    .disclaimer-box {{
        background-color: #2c1a1a;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #dc3545;
        color: #ffbaba;
        font-size: 0.85rem;
        margin-top: 30px;
    }}
    
    /* Premium Prompt text styling for the results section */
    .premium-prompt-text {{
        color: #ffffff;
        font-size: 0.95rem;
        margin-bottom: 8px;
    }}
    </style>
    """, unsafe_allow_html=True)

# --- 4. TOP NAVIGATION (Sign In & Go Premium) ---
nav_left, nav_right = st.columns([0.5, 0.5])
with nav_left:
    st.markdown("<h1 class='medical-title'>⚕️ CuraBot</h1>", unsafe_allow_html=True)
    st.markdown("<p class='medical-subtitle'>Intelligent Evidence-Based Health Insights</p>", unsafe_allow_html=True)

with nav_right:
    c1, c2 = st.columns(2)
    # AUTHENTICATION PANEL: Handles Login and Account Creation
    with c1:
        if not st.session_state.is_logged_in:
            with st.popover("👤 Sign In", use_container_width=True, disabled=st.session_state.loading):
                tab1, tab2 = st.tabs(["Login", "Create Account"])
                with tab1:
                    lemail = st.text_input("Email", key="l_email")
                    lpass = st.text_input("Password", type="password", key="l_pass")
                    if st.button("Log In", use_container_width=True):
                        if is_valid_email(lemail) and lpass:
                            st.session_state.is_logged_in = True
                            st.rerun()
                with tab2:
                    remail = st.text_input("Email", key="r_email")
                    rpass = st.text_input("Password", type="password", key="r_pass")
                    rconf = st.text_input("Confirm Password", type="password")
                    if st.button("Register", use_container_width=True):
                        if rpass == rconf and is_valid_email(remail) and len(rpass) >= 8:
                            st.session_state.is_logged_in = True
                            st.rerun()
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#28a745;'>✅ <b>Logged In</b></div>", unsafe_allow_html=True)

    # PREMIUM UPGRADE PANEL: Lists features and pricing
    with c2:
        if not st.session_state.is_premium:
            with st.popover("👑 Go Premium", use_container_width=True, disabled=st.session_state.loading):
                st.markdown("### 🔐 Premium Features")
                st.markdown('<div class="premium-feature">📄 <b>Download Reports:</b> Save PDF health reports.</div>', unsafe_allow_html=True)
                st.markdown('<div class="premium-feature">📊 <b>Comparison View:</b> Advanced condition analysis.</div>', unsafe_allow_html=True)
                st.markdown('<div class="premium-feature">🗂️ <b>History Tracking:</b> Secure long-term tracking.</div>', unsafe_allow_html=True)
                
                st.markdown("---")
                plan = st.radio("Choose Plan", ["Monthly - $9.99", "Yearly - $89.99 (Save 25%)"])
                
                if st.button("Upgrade to Premium", use_container_width=True):
                    if not st.session_state.is_logged_in:
                        st.warning("Please sign in to continue with Premium.")
                    else:
                        st.session_state.is_premium = True
                        st.success("Premium Activated!")
                        time.sleep(1)
                        st.rerun()
                st.caption("Maybe Later")
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#d4af37;'>⭐ <b>Premium Active</b></div>", unsafe_allow_html=True)

# --- 5. MAIN CONTENT: ASSESSMENT FORM ---
left_panel, right_panel = st.columns([2, 1])

with left_panel:
    # High-visibility Emergency Triage for user safety
    with st.expander("🚨 EMERGENCY TRIAGE", expanded=False):
        st.error("Seek immediate care for severe symptoms such as chest pain or breathing difficulties.")

    st.markdown("### 🔍 Symptom Assessment")
    user_input = st.text_area("Describe your symptoms:", height=180, disabled=st.session_state.loading, placeholder="E.g., Severe headache with light sensitivity...")

    # Demographic Inputs used for MeSH filtering and clinical context
    p1, p2, p3 = st.columns(3)
    with p1: age = st.selectbox("Age Group", ["Child (0-12)", "Teen (13-17)", "Adult (18-64)", "Senior (65+)"], index=2, disabled=st.session_state.loading)
    with p2: gender = st.selectbox("Gender", ["Male", "Female", "Other"], disabled=st.session_state.loading)
    with p3: severity = st.select_slider("Severity", options=["Mild", "Moderate", "Severe"], disabled=st.session_state.loading)

    # Action buttons: Analysis Trigger and Cancel
    b1, b2 = st.columns([0.7, 0.3])
    with b1: analyze_btn = st.button("🚀 Analyze Symptoms", type="primary", disabled=st.session_state.loading, use_container_width=True)
    with b2:
        if st.button("✖ Cancel", use_container_width=True):
            st.session_state.loading = False
            st.rerun()

# --- 6. ANALYSIS LOGIC (5 MIN TIMEOUT) ---
if analyze_btn and user_input:
    st.session_state.loading = True
    st.rerun()

if st.session_state.loading:
    st.markdown("---")
    prog = st.progress(0)
    st.info("🧬 **Searching Clinical Databases...** (Deep research can take 1-3 minutes)")
    try:
        payload = {"symptoms": f"{user_input} ({age}, {gender}, {severity})"}
        
        # INCREASED TIMEOUT: Changed from 300 to 600 seconds to give agents 
        # plenty of time to finish PubMed synthesis without the frontend cutting them off.
        response = requests.get("http://127.0.0.1:8000/analyze", params=payload, timeout=600)
        
        if response.status_code == 200:
            st.session_state.report_content = response.json().get("analysis")
            if st.session_state.is_logged_in:
                st.session_state.history.append(f"{datetime.now().strftime('%H:%M')} - {user_input[:15]}...")
            st.session_state.loading = False
            st.rerun()
        else:
            st.session_state.loading = False
            st.error(f"⚠️ **Server Error:** (Code {response.status_code}). The agents encountered an issue.")
            
    except requests.exceptions.ReadTimeout:
        # Specific handling for when the research exceeds 10 minutes
        st.session_state.loading = False
        st.error("⏱️ **Research Timeout:** The clinical research is taking an exceptional amount of time. Please try a more specific symptom description.")
    except Exception as e:
        st.session_state.loading = False
        st.error(f"❌ Connection Failed: {e}")

# --- 7. SIDEBAR & RESULTS DISPLAY ---
with right_panel:
    st.markdown("### 🔒 Privacy First")
    st.info("Symptoms remain anonymous. No personal data stored.")
    if not st.session_state.is_logged_in:
        st.markdown('<div style="background-color:#1e2129; padding:12px; border-radius:8px; border:1px dashed #4b9fff; font-size:0.85rem; color:#4b9fff;">💡 <b>Optional:</b> Sign in to save research history.</div>', unsafe_allow_html=True)
    if st.session_state.is_logged_in:
        st.markdown("### 🕒 Recent Activity")
        for h in st.session_state.history[-5:]: st.caption(f"📅 {h}")

# DISPLAY RESULTS AND GATED PREMIUM FEATURES
if st.session_state.report_content:
    st.markdown("---")
    res_col, pdf_col = st.columns([0.65, 0.35])
    
    with res_col: 
        st.markdown("### 📋 Clinical Findings")
    
    with pdf_col:
        # OPTION B: Clickable Upgrade Prompt with Text + Small Button
        if not st.session_state.is_premium:
            st.markdown('<p class="premium-prompt-text">👑 Want to save or share this report?</p>', unsafe_allow_html=True)
            # Full phrase "Upgrade to Premium" is the clickable action
            if st.button("✨ Upgrade to Premium", use_container_width=True):
                st.info("Please use the 'Go Premium' button in the top navigation bar to complete your upgrade.")
        else:
            # Full feature access for Premium users
            st.button("📄 Download PDF Report", use_container_width=True)

    # The AI-generated clinical summary
    st.markdown(f'<div class="result-card">{st.session_state.report_content}</div>', unsafe_allow_html=True)
    
    # Required Medical Disclaimer
    st.markdown("""
        <div class="disclaimer-box">
            <b>⚠️ MEDICAL DISCLAIMER:</b> Educational purposes only. This AI-synthesized report is not a substitute 
            for professional medical advice, diagnosis, or treatment. Always consult a physician.
        </div>
    """, unsafe_allow_html=True)
    