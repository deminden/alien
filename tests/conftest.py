from __future__ import annotations

import socket

import pytest


@pytest.fixture(autouse=True)
def _block_network(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail any test that tries to open a real network connection unless it is marked ``network``."""
    if request.node.get_closest_marker("network"):
        return

    original_connect = socket.socket.connect

    def guard(sock: socket.socket, address: object) -> None:
        if sock.family == getattr(socket, "AF_UNIX", None):
            original_connect(sock, address)  # process pools talk over local sockets
            return
        raise RuntimeError("Network access is disabled in the test suite; mark the test with @pytest.mark.network.")

    monkeypatch.setattr(socket.socket, "connect", guard)
