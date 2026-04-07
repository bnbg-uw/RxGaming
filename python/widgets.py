from __future__ import annotations

"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

widgets.py
Extensions and additions to PySide6 Qt widgets for the RxGaming tool.

This is a collection of generally unrelated widgets and widget helpers. They all extend Qt
in some way and act like Qt elements, so they are here together.
"""

from enum import Enum, auto
from typing import Callable

from PySide6.QtCore import Qt, QPoint
from PySide6.QtGui import QCloseEvent, QFontMetrics, QPaintEvent, QPainter, QPen
from PySide6.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QPushButton, QFileDialog, QSlider, QMainWindow


WindowCloseCallback = Callable[["QWindow"], None]
MainWindowCloseCallback = Callable[["QMainWindowRx"], None]


class QFileSelectionLineEdit(QWidget):
    """Line edit together with browse button to select files or directories."""

    class FileType(Enum):
        """Enumeration for indicating whether to open a file or a directory."""
        
        File = auto()
        """Open a file."""
        Directory = auto()
        """Open a directory."""
    
    def __init__(
        self,
        caption: str = "Browse...",
        file_type: QFileSelectionLineEdit.FileType = FileType.File,
        filter: str = "Any file (*.*)",
        new_file: bool = False,
    ) -> None:
        """Create a new ``QFileSelectionLineEdit``
        
        The file selection dialog window title is set by ``caption``. Pass in ``file_type`` to
        select whether to browse for files or directories. When browsing for files, ``filter``
        can be set to a standard filter text string (or list of strings) to limit possible choices.
        """
        
        super().__init__()
        self.caption = caption
        self.file_type = file_type
        self.filter = filter
        self.new_file = new_file
        
        self.main_layout = QHBoxLayout()
        self.line_edit = QLineEdit()
        self.button = QPushButton("Browse")
        self.main_layout.addWidget(self.line_edit)
        self.main_layout.addWidget(self.button)
        self.setLayout(self.main_layout)
        
        self.button.clicked.connect(self._on_browse)
    
    def _on_browse(self, checked: bool = False) -> None:
        # TODO: start from last place ended up
        if self.file_type is self.FileType.File:
            if self.new_file:
                file_path = QFileDialog.getSaveFileName(self, self.caption, "", self.filter)[0]
            else:
                file_path = QFileDialog.getOpenFileName(self, self.caption, "", self.filter)[0]
        else:
            file_path = QFileDialog.getExistingDirectory(self, self.caption, "")
        if file_path != "":
            self.line_edit.setText(file_path)
    
    def text(self) -> str:
        """Get text from the line edit."""
        
        return self.line_edit.text()
    
    def set_text(self, text: str) -> None:
        """Set text in the line edit."""
        self.line_edit.setText(text)


class QWindow(QWidget):
    """Extends QWidget for the purpose of hooking in to the window closure callback.
    
    To be used exactly like a QWidget that is a window (has no parent).
    """
    
    def __init__(self) -> None:
        super().__init__()
        self.on_closed: WindowCloseCallback | None = None
    
    def set_on_closed(self, callback: WindowCloseCallback) -> None:
        """Set callback for window close event.
        
        The ``callback`` should accept one argument, the ``QWindow`` that is closing.
        """
        
        self.on_closed = callback
    
    def closeEvent(self, event: QCloseEvent) -> None:
        """"""
        super().closeEvent(event)
        if self.on_closed is not None:
            self.on_closed(self)


class QMainWindowRx(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.on_closed: MainWindowCloseCallback | None = None

    def set_on_closed(self, callback: MainWindowCloseCallback) -> None:
        self.on_closed = callback

    def closeEvent(self, event: QCloseEvent) -> None:
        super().closeEvent(event)
        if self.on_closed is not None:
            self.on_closed(self)


class SliderWithValue(QSlider):
    def __init__(self, parent: QWidget | Qt.Orientation | None = None):
        super().__init__(parent)

        self.style_sheet = """
        QSlider::groove:vertical {
                background-color: #222;
                width: 30px;
        }
        QSlider::handle:vertical {
            border: 1px #438f99;
            border-style: outset;
            margin: -2px 0;
            width: 30px;
            height: 3px;
            background-color: #438f99;
        }
        QSlider::sub-page:vertical {
            background: #4B4B4B;
        }
        QSlider::groove:horizontal {
                background-color: #222;
                height: 30px;
        }
        QSlider::handle:horizontal {
            border: 1px #438f99;
            border-style: outset;
            margin: -2px 0;
            width: 3px;
            height: 30px;
            background-color: #438f99;
        }
        QSlider::sub-page:horizontal {
            background: #4B4B4B;
        }
        """

        # self.setStyleSheet(self.style_sheet)

    def paintEvent(self, event: QPaintEvent) -> None:
        super().paintEvent(event)

        curr_value = str(self.value())
        round_value = round(float(curr_value), 2)

        painter = QPainter(self)
        painter.setPen(QPen(Qt.GlobalColor.black))

        font_metrics = QFontMetrics(self.font())
        font_width = font_metrics.boundingRect(str(round_value)).width()
        rect = self.rect()
        if self.orientation() == Qt.Orientation.Horizontal:
            horizontal_x_pos = rect.width() - font_width - 5
            horizontal_y_pos = int(rect.height() * 0.75)

            painter.drawText(QPoint(horizontal_x_pos, horizontal_y_pos), str(round_value))

        elif self.orientation() == Qt.Orientation.Vertical:
            painter.drawText(
                QPoint(int(rect.width() / 2.0 - font_width / 2.0), rect.height() - 5),
                str(round_value),
            )

        painter.drawRect(rect)
