from fastapi import FastAPI, HTTPException
from crewai import Crew, Process
from backend.agents import classifier_agent, matcher_agent, advisor_agent
from backend.tasks import create_tasks
import time

# Initialize FastAPI app
app = FastAPI(title="CuraBot Backend")

@app.get("/analyze")
async def analyze(symptoms: str):
    max_retries = 2
    for attempt in range(max_retries):
        try:
            tasks = create_tasks(symptoms, classifier_agent, matcher_agent, advisor_agent)
            
            crew = Crew(
                agents=[classifier_agent, matcher_agent, advisor_agent],
                tasks=tasks,
                process=Process.sequential,
                max_rpm=3, # Lowered to 3 to stay safely within free-tier limits
                verbose=True,
                cache=True
            )
            
            result = crew.kickoff()
            final_output = getattr(result, 'raw', str(result))
            return {"analysis": final_output}
            
        except Exception as e:
            if "rate_limit" in str(e).lower() and attempt < max_retries - 1:
                time.sleep(15) # Cooldown before retry
                continue
            raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    # Use 'timeout_keep_alive' to ensure the Uvicorn worker doesn't kill the process mid-research
    uvicorn.run(app, host="0.0.0.0", port=8000, timeout_keep_alive=600)