// OCSR Service using DECIMER
import api from "./api";
import type { OCSRResult } from "../types/api";

const OCSR_URL = "/ocsr";

export const processImage = async (
  image: File | Blob,
  reference: string | null = null,
  handDrawn = false
): Promise<OCSRResult> => {
  const formData = new FormData();
  formData.append("file", image);
  formData.append("hand_drawn", String(handDrawn));
  if (reference) {
    formData.append("reference", reference);
  }

  try {
    const response = await api.post<OCSRResult>(`${OCSR_URL}/process-upload`, formData, {
      headers: {
        "Content-Type": "multipart/form-data",
      },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`OCSR processing failed: ${message}`);
  }
};

export const processImageUrl = async (
  imageUrl: string,
  reference: string | null = null,
  handDrawn = false
): Promise<OCSRResult> => {
  try {
    const response = await api.post<OCSRResult>(`${OCSR_URL}/process`, {
      path: imageUrl,
      reference: reference || undefined,
      hand_drawn: handDrawn,
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`OCSR processing failed: ${message}`);
  }
};

const ocsrService = {
  processImage,
  processImageUrl,
};

export default ocsrService;
