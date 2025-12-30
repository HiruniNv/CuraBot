from fastapi import FastAPI, HTTPException, Request, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from crewai import Crew, Process
from backend.agents import classifier_agent, matcher_agent, advisor_agent
from backend.tasks import create_tasks
import hashlib
import os
import time

# ============================================================================
# CONFIGURATION
# ============================================================================
PAYHERE_MERCHANT_ID = os.getenv("PAYHERE_MERCHANT_ID", "1226123") 
PAYHERE_MERCHANT_SECRET = os.getenv("PAYHERE_MERCHANT_SECRET", "") 

app = FastAPI(title="CuraBot Backend API", version="2.0.4")

# CORS Setup
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # In production, replace with ["http://localhost:8501"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class PaymentRequest(BaseModel):
    plan: str
    user_email: str
    user_name: str = "CuraBot User"

# ============================================================================
# AI ANALYSIS ENDPOINT
# ============================================================================
@app.get("/analyze")
async def analyze(symptoms: str = Query(..., min_length=3)):
    try:
        # Create tasks based on the symptoms
        tasks = create_tasks(symptoms, classifier_agent, matcher_agent, advisor_agent)
        
        # Initialize the Crew
        # max_rpm=2 ensures the whole crew stays within Groq Free Tier limits
        crew = Crew(
            agents=[classifier_agent, matcher_agent, advisor_agent],
            tasks=tasks,
            process=Process.sequential,
            verbose=True,
            max_rpm=2 
        )
        
        # Kickoff the research process
        result = crew.kickoff()
        
        # Robust extraction of string output
        final_output = result.raw if hasattr(result, 'raw') else str(result)

        return {"analysis": final_output}
        
    except Exception as e:
        error_str = str(e).lower()
        print(f"❌ Analysis failed: {error_str}")
        
        if "rate_limit_exceeded" in error_str or "429" in error_str:
            raise HTTPException(
                status_code=429, 
                detail="CuraBot's research engine is currently at capacity. Please wait 60 seconds."
            )
            
        raise HTTPException(status_code=500, detail="AI research failed. Please check backend logs.")

# ============================================================================
# PAYMENT ENDPOINT (PAYHERE)
# ============================================================================
@app.post("/create-payment")
async def create_payment(request: PaymentRequest):
    try:
        # Amount logic: Monthly vs Yearly
        amount_val = 999.00 if "Monthly" in request.plan else 8999.00
        formatted_amount = "{:.2f}".format(amount_val) 
        
        order_id = f"CURA{int(time.time())}"
        currency = "LKR"
        
        # --- PAYHERE HASH GENERATION ---
        # Uppercase(Md5( MerchantID + OrderID + Amount + Currency + Uppercase(Md5(Secret)) ))
        secret_md5 = hashlib.md5(PAYHERE_MERCHANT_SECRET.encode()).hexdigest().upper()
        hash_payload = f"{PAYHERE_MERCHANT_ID}{order_id}{formatted_amount}{currency}{secret_md5}"
        merchant_hash = hashlib.md5(hash_payload.encode()).hexdigest().upper()
        
        return {
            "merchant_id": PAYHERE_MERCHANT_ID,
            "return_url": "http://localhost:8501/?payment=success",
            "cancel_url": "http://localhost:8501/?payment=cancel",
            "notify_url": "https://your-backend-live-url.com/payment-notify",
            "order_id": order_id,
            "items": f"CuraBot Premium: {request.plan}",
            "currency": currency,
            "amount": formatted_amount,
            "first_name": request.user_name.split()[0] if request.user_name else "Guest",
            "last_name": "User",
            "email": request.user_email,
            "hash": merchant_hash,
            # Sandbox toggle: Set to True for testing
            "sandbox": True 
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Payment error: {str(e)}")

# ============================================================================
# SERVER START
# ============================================================================
if __name__ == "__main__":
    import uvicorn
    # CRITICAL: Increased timeout to 300s to prevent 'Gateway Timeout' during AI generation
    uvicorn.run(app, host="0.0.0.0", port=8000, timeout_keep_alive=300)