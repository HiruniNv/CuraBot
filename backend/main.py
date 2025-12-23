from fastapi import FastAPI, HTTPException
from crewai import Crew, Process
from .agents import classifier_agent, matcher_agent, advisor_agent
from .tasks import create_tasks

app = FastAPI(title="CuraBot Backend")

@app.get("/analyze")
async def analyze(symptoms: str):
    try:
        tasks = create_tasks(symptoms, classifier_agent, matcher_agent, advisor_agent)
        crew = Crew(
            agents=[classifier_agent, matcher_agent, advisor_agent],
            tasks=tasks,
            process=Process.sequential
        )
        result = crew.kickoff()
        return {"result": str(result)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)