#app.py

import streamlit as st
import streamlit.components.v1 as components
import requests
import pandas as pd
from fpdf import FPDF
from datetime import datetime
import time
import re

# --- 1. PAGE CONFIGURATION & SESSION STATE ---
st.set_page_config(
    page_title="CuraBot | AI Medical Research", 
    page_icon="⚕️", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# PAYHERE PAYMENT SUCCESS HANDLER
# ============================================================================
query_params = st.query_params
if query_params.get("payment") == "success":
    st.session_state.is_premium = True
    st.success("🎉 Payment Successful! Welcome to CuraBot Premium.")
    st.query_params.clear()
elif query_params.get("payment") == "cancel":
    st.warning("⚠️ Payment was cancelled. You can try again anytime.")
    st.query_params.clear()

# Initialize Session States
states = {
    'is_logged_in': False,
    'is_premium': False,
    'history': [],
    'loading': False,
    'report_content': None,
    'user_email': None
}
for key, val in states.items():
    if key not in st.session_state:
        st.session_state[key] = val

# --- 2. HELPER FUNCTIONS ---
def is_valid_email(email):
    return re.match(r"[^@]+@[^@]+\.[^@]+", email)

def generate_pdf(report_text, symptoms, demographics):
    """Creates a professional PDF report with clinical findings."""
    pdf = FPDF()
    pdf.add_page()
    
    # Header
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, "CuraBot Clinical Research Report", ln=True, align='C')
    pdf.set_font("Arial", 'I', 10)
    pdf.cell(200, 10, f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M')}", ln=True, align='C')
    pdf.ln(10)
    
    # Patient Context
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Case Summary:", ln=True)
    pdf.set_font("Arial", '', 11)
    pdf.multi_cell(0, 10, f"Symptoms: {symptoms}\nContext: {demographics}")
    pdf.ln(5)
    
    # Findings
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, "Clinical Synthesis:", ln=True)
    pdf.set_font("Arial", '', 11)
    clean_text = report_text.replace("**", "")
    pdf.multi_cell(0, 8, clean_text)
    
    # Disclaimer
    pdf.ln(10)
    pdf.set_font("Arial", 'I', 8)
    pdf.multi_cell(0, 5, "DISCLAIMER: This report is AI-generated for educational purposes and does not constitute medical advice.")
    
    return pdf.output(dest='S').encode('latin-1', 'replace')

# ============================================================================
# PREMIUM DIALOG FUNCTION (Updated for centered modal)
# ============================================================================
@st.dialog("👑 Upgrade to CuraBot Premium", width="large")
def show_premium_dialog():
    """Display premium upgrade options in a centered dialog modal"""
    
    if not st.session_state.is_logged_in:
        st.warning("⚠️ Please Sign In first to upgrade.")
        return
    
    st.markdown("### Unlock Professional Health Intelligence")
    st.markdown("Get access to advanced features designed for comprehensive health management.")
    
    # Feature showcase with visual styling
    st.markdown('<div class="premium-feature">📄 <b>Download PDF Reports:</b> Save and share professional clinical reports</div>', unsafe_allow_html=True)
    st.markdown('<div class="premium-feature">📊 <b>Advanced Analytics:</b> Detailed condition comparison and risk assessment</div>', unsafe_allow_html=True)
    st.markdown('<div class="premium-feature">🔍 <b>Extended Research:</b> Access to full PubMed literature analysis</div>', unsafe_allow_html=True)
    st.markdown('<div class="premium-feature">📈 <b>Health History:</b> Track and compare symptoms over time</div>', unsafe_allow_html=True)
    
    st.markdown("---")
    
    # Plan selection
    plan = st.radio(
        "Choose Your Plan", 
        ["Monthly - Rs. 999", "Yearly - Rs. 8,999 (Save 25%)"], 
        key="dialog_plan",
        help="Annual plan saves you Rs. 3,000 per year!"
    )
    
    if st.button("🚀 Proceed to Secure Checkout", use_container_width=True, type="primary"):
        with st.spinner("Initializing secure payment gateway..."):
            try:
                res = requests.post(
                    "http://127.0.0.1:8000/create-payment",
                    json={
                        "plan": plan,
                        "user_email": st.session_state.user_email,
                        "user_name": "CuraBot User"
                    },
                    timeout=30
                )
                
                if res.status_code == 200:
                    p = res.json()
                    
                    # Auto-submit form to PayHere with target="_top" to break iframe
                    form_html = f"""
                    <html>
                    <head><title>Redirecting to Payment...</title></head>
                    <body>
                        <form method="post" action="https://sandbox.payhere.lk/pay/checkout" id="payform" target="_top">
                            <input type="hidden" name="merchant_id" value="{p['merchant_id']}">
                            <input type="hidden" name="return_url" value="{p['return_url']}">
                            <input type="hidden" name="cancel_url" value="{p['cancel_url']}">
                            <input type="hidden" name="notify_url" value="{p['notify_url']}">
                            <input type="hidden" name="order_id" value="{p['order_id']}">
                            <input type="hidden" name="items" value="{p['items']}">
                            <input type="hidden" name="currency" value="LKR">
                            <input type="hidden" name="amount" value="{p['amount']}">
                            <input type="hidden" name="first_name" value="{p['first_name']}">
                            <input type="hidden" name="last_name" value="User">
                            <input type="hidden" name="email" value="{p['email']}">
                            <input type="hidden" name="hash" value="{p['hash']}">
                        </form>
                        <script type="text/javascript">
                            window.onload = function() {{
                                document.getElementById('payform').submit();
                            }};
                        </script>
                    </body>
                    </html>
                    """
                    components.html(form_html, height=0)
                    st.success("✅ Payment gateway initialized!")
                    st.info("🔄 Redirecting to PayHere secure checkout...")
                else:
                    st.error(f"❌ Payment initialization failed: {res.text}")
                    
            except requests.exceptions.Timeout:
                st.error("⚠️ Connection Timeout: Payment server is taking too long to respond. Please try again.")
            except Exception as e:
                st.error(f"❌ Connection Error: {str(e)}")

# --- 3. CUSTOM CSS STYLING ---
st.markdown("""
    <style>
    /* Dark background */
    .main { background-color: #0e1117; }
    
    /* Branding */
    .medical-title { font-size: 2.8rem; font-weight: 800; color: #ffffff; margin-bottom: 0px; }
    .medical-subtitle { color: #4b9fff; font-size: 1.1rem; font-weight: 500; margin-bottom: 25px; }
    
    /* Primary Action Button */
    div.stButton > button:first-child[kind="primary"] {
        background-color: #28a745 !important;
        color: white !important;
        border: none !important;
        height: 3.5em !important;
        width: 100% !important;
        font-weight: bold !important;
        font-size: 1.1rem !important;
    }
    
    /* Premium Features in Dialog */
    .premium-feature {
        padding: 12px 15px;
        border-radius: 8px;
        background: #1e2129;
        margin-bottom: 10px;
        border-left: 4px solid #d4af37;
        font-size: 0.95rem;
    }

    /* Results Card */
    .result-card {
        background-color: #1e2129;
        padding: 25px;
        border-radius: 15px;
        border-left: 5px solid #4b9fff;
        margin-top: 20px;
        line-height: 1.6;
    }

    /* Disclaimer Box */
    .disclaimer-box {
        background-color: #2c1a1a;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #dc3545;
        color: #ffbaba;
        font-size: 0.85rem;
        margin-top: 30px;
    }
    
    /* Premium Prompt Text */
    .premium-prompt-text {
        color: #ffffff;
        font-size: 0.95rem;
        margin-bottom: 8px;
        font-weight: 500;
    }
    
    /* Dialog Styling Override */
    div[data-testid="stDialog"] {
        background-color: #1e2129 !important;
    }
    </style>
    """, unsafe_allow_html=True)

# --- 4. TOP NAVIGATION (Updated with Dialog Trigger) ---
nav_left, nav_right = st.columns([0.5, 0.5])

with nav_left:
    st.markdown("<h1 class='medical-title'>⚕️ CuraBot</h1>", unsafe_allow_html=True)
    st.markdown("<p class='medical-subtitle'>Intelligent Evidence-Based Health Insights</p>", unsafe_allow_html=True)

with nav_right:
    c1, c2 = st.columns(2)
    
    # AUTHENTICATION PANEL
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
                            st.session_state.user_email = lemail
                            st.rerun()
                        else:
                            st.error("Invalid email or password")
                
                with tab2:
                    remail = st.text_input("Email", key="r_email")
                    rpass = st.text_input("Password", type="password", key="r_pass")
                    rconf = st.text_input("Confirm Password", type="password")
                    if st.button("Register", use_container_width=True):
                        if rpass == rconf and is_valid_email(remail) and len(rpass) >= 8:
                            st.session_state.is_logged_in = True
                            st.session_state.user_email = remail
                            st.success("Account created successfully!")
                            st.rerun()
                        else:
                            st.error("Please ensure passwords match and are at least 8 characters")
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#28a745;'>✅ <b>Logged In</b></div>", unsafe_allow_html=True)

    # PREMIUM UPGRADE BUTTON (Triggers Dialog)
    with c2:
        if not st.session_state.is_premium:
            if st.button("👑 Go Premium", use_container_width=True, disabled=st.session_state.loading):
                show_premium_dialog()
        else:
            st.markdown("<div style='text-align:center; padding:10px; color:#d4af37;'>⭐ <b>Premium Active</b></div>", unsafe_allow_html=True)

# --- 5. MAIN CONTENT: ASSESSMENT FORM ---
left_panel, right_panel = st.columns([2, 1])

with left_panel:
    # Emergency Triage
    with st.expander("🚨 EMERGENCY TRIAGE", expanded=False):
        st.error("**Seek immediate care for:**")
        st.markdown("""
        - Severe chest pain or pressure
        - Difficulty breathing or shortness of breath
        - Sudden severe headache
        - Loss of consciousness
        - Severe allergic reactions
        - Uncontrolled bleeding
        """)

    st.markdown("### 🔍 Symptom Assessment")
    user_input = st.text_area(
        "Describe your symptoms in detail:",
        height=180,
        disabled=st.session_state.loading,
        placeholder="E.g., I've been experiencing severe headaches with sensitivity to light for the past 3 days. The pain is throbbing and gets worse with movement..."
    )

    # Demographic Inputs
    p1, p2, p3 = st.columns(3)
    with p1:
        age = st.selectbox(
            "Age Group",
            ["Child (0-12)", "Teen (13-17)", "Adult (18-64)", "Senior (65+)"],
            index=2,
            disabled=st.session_state.loading
        )
    with p2:
        gender = st.selectbox(
            "Gender",
            ["Male", "Female", "Other"],
            disabled=st.session_state.loading
        )
    with p3:
        severity = st.select_slider(
            "Severity",
            options=["Mild", "Moderate", "Severe"],
            disabled=st.session_state.loading
        )

    # Action buttons
    b1, b2 = st.columns([0.7, 0.3])
    with b1:
        analyze_btn = st.button(
            "🚀 Analyze Symptoms",
            type="primary",
            disabled=st.session_state.loading or not user_input,
            use_container_width=True
        )
    with b2:
        if st.button("✖ Cancel", use_container_width=True, disabled=not st.session_state.loading):
            st.session_state.loading = False
            st.rerun()

# --- 6. ANALYSIS LOGIC ---
if analyze_btn and user_input:
    st.session_state.loading = True
    st.rerun()

if st.session_state.loading:
    st.markdown("---")
    prog = st.progress(0)
    status_text = st.empty()
    
    status_text.info("🧬 **Searching Clinical Databases...** (This may take 1-3 minutes for comprehensive analysis)")
    
    try:
        payload = {"symptoms": f"{user_input} ({age}, {gender}, {severity})"}
        
        # 10-minute timeout for deep research
        response = requests.get(
            "http://127.0.0.1:8000/analyze",
            params=payload,
            timeout=600,
            stream=True
        )
        
        # Simulate progress for better UX
        for i in range(100):
            prog.progress(i + 1)
            time.sleep(0.02)
        
        if response.status_code == 200:
            st.session_state.report_content = response.json().get("analysis")
            
            # Log activity for signed-in users
            if st.session_state.is_logged_in:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
                entry = f"{timestamp} - {user_input[:50]}{'...' if len(user_input) > 50 else ''}"
                st.session_state.history.append(entry)
            
            st.session_state.loading = False
            status_text.success("✅ Analysis complete!")
            time.sleep(1)
            st.rerun()
        else:
            st.session_state.loading = False
            st.error(f"⚠️ **Server Error** (Code {response.status_code}): The analysis agents encountered an issue. Please try again.")
            
    except requests.exceptions.ReadTimeout:
        st.session_state.loading = False
        st.error("⏱️ **Research Timeout:** The clinical research took longer than expected. Try using more specific symptom descriptions.")
    except requests.exceptions.ConnectionError:
        st.session_state.loading = False
        st.error("❌ **Connection Failed:** Cannot reach the analysis server. Please ensure the backend is running.")
    except Exception as e:
        st.session_state.loading = False
        st.error(f"❌ **Unexpected Error:** {str(e)}")

# --- 7. SIDEBAR & ACTIVITY PANEL ---
with right_panel:
    st.markdown("### 🔒 Privacy First")
    st.info("Your symptoms remain completely anonymous. We don't store any personal health information.")
    
    if not st.session_state.is_logged_in:
        st.markdown("""
        <div style="background-color:#1e2129; padding:12px; border-radius:8px; 
                    border:1px dashed #4b9fff; font-size:0.85rem; color:#4b9fff;">
            💡 <b>Tip:</b> Sign in to save your research history and track symptoms over time.
        </div>
        """, unsafe_allow_html=True)
    
    if st.session_state.is_logged_in and st.session_state.history:
        st.markdown("### 🕒 Recent Activity")
        for h in st.session_state.history[-5:]:
            st.caption(f"📅 {h}")

# ============================================================================
# RESULTS DISPLAY WITH PREMIUM GATING
# ============================================================================
if st.session_state.report_content:
    st.markdown("---")
    res_col, pdf_col = st.columns([0.65, 0.35])
    
    with res_col:
        st.markdown("### 📋 Clinical Findings")
    
    with pdf_col:
        if not st.session_state.is_premium:
            # Show premium upgrade prompt
            st.markdown('<p class="premium-prompt-text">👑 Want to save this report?</p>', unsafe_allow_html=True)
            if st.button("✨ Upgrade to Premium", use_container_width=True, key="upgrade_from_results"):
                show_premium_dialog()
        else:
            # Enable PDF download for premium users
            demo_string = f"{age}, {gender}, {severity} severity"
            pdf_data = generate_pdf(st.session_state.report_content, user_input, demo_string)
            
            st.download_button(
                label="📄 Download PDF Report",
                data=pdf_data,
                file_name=f"CuraBot_Report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                mime="application/pdf",
                use_container_width=True
            )

    # Display AI-generated clinical report
    st.markdown(
        f'<div class="result-card">{st.session_state.report_content}</div>',
        unsafe_allow_html=True
    )
    
    # Medical Disclaimer
    st.markdown("""
        <div class="disclaimer-box">
            <b>⚠️ MEDICAL DISCLAIMER:</b> This analysis is provided for educational and informational purposes only.
            The AI-generated content does not constitute professional medical advice, diagnosis, or treatment.
            Always seek the advice of qualified health providers with any questions regarding medical conditions.
            Never disregard professional medical advice or delay seeking it because of information provided by this tool.
        </div>
    """, unsafe_allow_html=True)

# --- 8. FOOTER ---
st.markdown("---")
footer_cols = st.columns([1, 1, 1])
with footer_cols[0]:
    st.caption("🏥 Powered by Advanced AI Research")
with footer_cols[1]:
    st.caption("🔬 Evidence-Based Clinical Intelligence")
with footer_cols[2]:
    st.caption(f"📅 {datetime.now().strftime('%Y')} CuraBot")