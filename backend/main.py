from fastapi import FastAPI, HTTPException
from crewai import Crew, Process
from backend.agents import classifier_agent, matcher_agent, advisor_agent
from backend.tasks import create_tasks
import time

app = FastAPI(title="CuraBot Backend")

@app.get("/analyze")
async def analyze(symptoms: str):
    max_retries = 2
    for attempt in range(max_retries):
        try:
            # Generate tasks with current symptoms
            tasks = create_tasks(symptoms, classifier_agent, matcher_agent, advisor_agent)
            
            # Configure the Crew
            crew = Crew(
                agents=[classifier_agent, matcher_agent, advisor_agent],
                tasks=tasks,
                process=Process.sequential,
                max_rpm=3 # Keeps  API keys safe from rate limits
            )
            
            # Execute workflow
            result = crew.kickoff()
            
            # CrewAI 0.x returns a string or object; cast to string for the frontend
            return {"analysis": str(result)}
            
        except Exception as e:
            if "rate_limit" in str(e).lower() and attempt < max_retries - 1:
                time.sleep(20) # Wait for rate limit reset
                continue
            raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
    