from pydantic import BaseModel


class Msg(BaseModel):
    """
    A Pydantic model class representing a message.

    Attributes:
        msg (str): The message content.
    """

    msg: str
