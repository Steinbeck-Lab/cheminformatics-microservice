// OCSR Service using DECIMER
import { failed } from "@/lib/error-messages";
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
    throw failed("OCSR processing failed", error);
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
    throw failed("OCSR processing failed", error);
  }
};

const ocsrService = {
  processImage,
  processImageUrl,
};

export default ocsrService;
