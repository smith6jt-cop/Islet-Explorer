"""
AI Assistant component for Islet Explorer Panel App.
Integrates with UF Navigator Toolkit API for chat functionality.
"""

import os
from typing import Optional, Dict, Any, Callable
import httpx
import panel as pn

# Navigator API configuration
NAVIGATOR_API_URL = os.environ.get(
    "NAVIGATOR_API_URL",
    "https://api.rc.ufl.edu/ai/ufos/v1/chat/completions"
)
NAVIGATOR_API_KEY = os.environ.get("NAVIGATOR_API_KEY", "")

# Default models
MODELS = {
    "Navigator Fast (20B)": "gpt-oss-20b",
    "Navigator Large (210B)": "gpt-oss-210b"
}

# System prompts for different contexts
SYSTEM_PROMPTS = {
    "islet_analysis": """You are an expert in pancreatic islet biology and spatial omics analysis.
You help researchers interpret data from multiplex imaging of human pancreatic tissue.
You understand:
- Islet cell types (alpha, beta, delta cells) and their markers (glucagon/GCG, insulin/INS, somatostatin/SST)
- Disease progression from normal (ND) through autoantibody positive (Aab+) to Type 1 Diabetes (T1D)
- Spatial relationships within islets (core vs band regions)
- Statistical analysis of cell populations and marker expression
- Islet size distributions and their relationship to function

Provide clear, scientifically accurate interpretations. When discussing statistics, explain both the methods and biological significance.""",

    "statistical": """You are a biostatistician specializing in spatial omics data analysis.
You help researchers understand:
- ANOVA and post-hoc comparisons for group differences
- Effect sizes and clinical significance
- Multiple testing correction (Benjamini-Hochberg)
- Distribution analysis and normality testing
- Outlier detection and handling
- Power analysis and sample size considerations

Provide practical guidance on statistical interpretation in the context of biomedical research.""",

    "troubleshooting": """You are a technical support assistant for the Islet Explorer application.
You help users with:
- Data loading and formatting issues
- Understanding plot controls and options
- Interpreting error messages
- Optimizing analysis workflows
- Export and reporting features

Be concise and provide step-by-step solutions when applicable."""
}


class NavigatorClient:
    """Client for UF Navigator AI API."""

    def __init__(
        self,
        api_url: str = NAVIGATOR_API_URL,
        api_key: str = NAVIGATOR_API_KEY
    ):
        self.api_url = api_url
        self.api_key = api_key
        self.total_tokens = 0

    async def chat_completion(
        self,
        messages: list,
        model: str = "gpt-oss-20b",
        temperature: float = 0.7,
        max_tokens: int = 1024
    ) -> Dict[str, Any]:
        """Send chat completion request to Navigator API."""

        if not self.api_key:
            return {
                "error": "API key not configured. Set NAVIGATOR_API_KEY environment variable.",
                "content": None
            }

        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }

        payload = {
            "model": model,
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens
        }

        try:
            async with httpx.AsyncClient(timeout=60.0) as client:
                response = await client.post(
                    self.api_url,
                    headers=headers,
                    json=payload
                )
                response.raise_for_status()
                result = response.json()

                # Track token usage
                if "usage" in result:
                    self.total_tokens += result["usage"].get("total_tokens", 0)

                return {
                    "content": result["choices"][0]["message"]["content"],
                    "usage": result.get("usage", {}),
                    "error": None
                }

        except httpx.TimeoutException:
            return {"error": "Request timed out. Please try again.", "content": None}
        except httpx.HTTPStatusError as e:
            return {"error": f"API error: {e.response.status_code}", "content": None}
        except Exception as e:
            return {"error": f"Error: {str(e)}", "content": None}


class AIAssistantPanel:
    """Panel-based AI Assistant chat interface."""

    def __init__(
        self,
        context: str = "islet_analysis",
        on_plot_context: Optional[Callable] = None
    ):
        self.client = NavigatorClient()
        self.context = context
        self.on_plot_context = on_plot_context
        self.chat_history = []

        # Create widgets
        self.model_selector = pn.widgets.Select(
            name="Model",
            options=list(MODELS.keys()),
            value="Navigator Fast (20B)",
            width=200
        )

        self.context_selector = pn.widgets.Select(
            name="Context",
            options=list(SYSTEM_PROMPTS.keys()),
            value=context,
            width=200
        )

        self.token_indicator = pn.pane.Markdown(
            "Tokens used: 0",
            styles={"font-size": "12px", "color": "#666"}
        )

        self.clear_button = pn.widgets.Button(
            name="Clear Chat",
            button_type="warning",
            width=100
        )
        self.clear_button.on_click(self._clear_chat)

        # Create chat interface
        self.chat_interface = pn.chat.ChatInterface(
            callback=self._chat_callback,
            user="You",
            avatar="🧑‍🔬",
            show_rerun=False,
            show_undo=False,
            sizing_mode="stretch_width"
        )

        # Add system message
        self.chat_interface.send(
            "Hello! I'm your AI assistant for islet analysis. "
            "Ask me about your data, statistical interpretations, or how to use the application.",
            user="Assistant",
            avatar="🤖",
            respond=False
        )

    async def _chat_callback(
        self,
        contents: str,
        user: str,
        instance: pn.chat.ChatInterface
    ):
        """Handle chat messages."""

        # Get system prompt
        system_prompt = SYSTEM_PROMPTS.get(
            self.context_selector.value,
            SYSTEM_PROMPTS["islet_analysis"]
        )

        # Add plot context if available
        if self.on_plot_context:
            plot_context = self.on_plot_context()
            if plot_context:
                system_prompt += f"\n\nCurrent analysis context:\n{plot_context}"

        # Build messages
        messages = [{"role": "system", "content": system_prompt}]

        # Add chat history
        for msg in self.chat_history[-10:]:  # Keep last 10 exchanges
            messages.append(msg)

        # Add current message
        messages.append({"role": "user", "content": contents})

        # Get model
        model = MODELS.get(self.model_selector.value, "gpt-oss-20b")

        # Make API call
        yield "Thinking..."

        result = await self.client.chat_completion(messages, model=model)

        if result["error"]:
            yield f"**Error:** {result['error']}"
            return

        response_content = result["content"]

        # Update chat history
        self.chat_history.append({"role": "user", "content": contents})
        self.chat_history.append({"role": "assistant", "content": response_content})

        # Update token counter
        self.token_indicator.object = f"Tokens used: {self.client.total_tokens:,}"

        yield response_content

    def _clear_chat(self, event):
        """Clear chat history."""
        self.chat_history.clear()
        self.chat_interface.clear()
        self.chat_interface.send(
            "Chat cleared. How can I help you?",
            user="Assistant",
            avatar="🤖",
            respond=False
        )

    def get_panel(self) -> pn.Column:
        """Return the complete chat panel."""
        header = pn.Row(
            pn.pane.Markdown("## AI Assistant", styles={"margin": "0"}),
            pn.Spacer(),
            self.clear_button,
            align="center"
        )

        controls = pn.Row(
            self.model_selector,
            self.context_selector,
            pn.Spacer(),
            self.token_indicator,
            sizing_mode="stretch_width"
        )

        return pn.Column(
            header,
            controls,
            self.chat_interface,
            sizing_mode="stretch_both",
            styles={
                "background": "#f8f9fa",
                "border-radius": "8px",
                "padding": "10px"
            }
        )


def create_ai_assistant(
    context: str = "islet_analysis",
    on_plot_context: Optional[Callable] = None
) -> pn.Column:
    """Create and return AI assistant panel."""
    assistant = AIAssistantPanel(context=context, on_plot_context=on_plot_context)
    return assistant.get_panel()
