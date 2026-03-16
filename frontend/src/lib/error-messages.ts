type DomainMessages = {
  default: string;
  parse?: string;
  timeout?: string;
};

const ERROR_MESSAGES: Record<string, DomainMessages> = {
  chem: {
    default: "Could not analyze this molecule. Please check the input and try again.",
    parse: "The molecule could not be parsed. Verify the SMILES or structure input.",
    timeout: "Analysis is taking too long. Try a simpler molecule or try again later.",
  },
  convert: {
    default: "Could not convert this structure. Please verify the input format.",
    parse: "The input could not be parsed for conversion. Check the format.",
    timeout: "Conversion timed out. Try a smaller molecule or try again later.",
  },
  depict: {
    default: "Could not generate the depiction. Please check the input molecule.",
    parse: "The molecule could not be parsed for depiction.",
    timeout: "Depiction generation timed out. Try a simpler structure.",
  },
  tools: {
    default: "The tool encountered an error. Please try again.",
    parse: "The input could not be processed by this tool.",
    timeout: "Processing timed out. Try again with simpler input.",
  },
  ocsr: {
    default: "Could not recognize the chemical structure from the image.",
    parse: "The image could not be processed. Try a clearer image.",
    timeout: "Image recognition timed out. Try a smaller or clearer image.",
  },
  network: {
    default: "Service temporarily unavailable. Please try again in a moment.",
    timeout: "The request timed out. Please check your connection and try again.",
  },
};

function isTimeoutError(error: unknown): boolean {
  if (error instanceof Error) {
    return error.message.includes("timeout") || error.message.includes("ECONNABORTED");
  }
  if (typeof error === "object" && error !== null && "code" in error) {
    return (error as { code: string }).code === "ECONNABORTED";
  }
  return false;
}

function isParseError(error: unknown): boolean {
  if (error instanceof Error) {
    return (
      error.message.includes("parse") ||
      error.message.includes("invalid") ||
      error.message.includes("422")
    );
  }
  if (typeof error === "object" && error !== null && "response" in error) {
    const resp = (error as { response?: { status?: number } }).response;
    return resp?.status === 422;
  }
  return false;
}

export function getErrorMessage(domain: string, error: unknown): string {
  const messages = ERROR_MESSAGES[domain] ?? ERROR_MESSAGES.network;

  if (import.meta.env.DEV) {
    console.error(`[${domain}] Error:`, error);
  }

  if (isTimeoutError(error)) {
    return messages.timeout ?? messages.default;
  }

  if (isParseError(error)) {
    return messages.parse ?? messages.default;
  }

  return messages.default;
}
